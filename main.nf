#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-

==============================
IPAW: HiRIEF II varDB pipeline
==============================
@Authors
Jorrit Boekel @glormph
Yafeng Zhu @yafeng

https://github.com/lehtiolab/proteogenomics-analysis-workflow
*/

nf_required_version = '19.04.0'
if( ! nextflow.version.matches(">= ${nf_required_version}") ){
  println("Nextflow version too old, ${nf_required_version} required")
  exit(1)
}


/* SET DEFAULT PARAMS */
params.isobaric = false
params.activation = 'hcd'
params.bamfiles = false
params.mods = false
params.novheaders = false
params.varheaders = false
params.saavheader = false
params.noclassfdr = false
params.dbsnp = false
params.cosmic = false
params.pisepdb = false
params.mzmldef = false
params.denoms = false
params.normalpsms = false
params.annovar_dir = false
params.bigwigs = false
params.splitchar = false
params.quantlookup = false

mods = file(params.mods)
knownproteins = file(params.knownproteins)
blastdb = file(params.blastdb)
gtffile = file(params.gtf)
snpfa = file(params.snpfa)
dbsnp = params.dbsnp ? file(params.dbsnp) : false
cosmic = params.cosmic ? file(params.cosmic) : false
genomefa = file(params.genome)
tdb = file(params.tdb)
normalpsms = params.normalpsms ? file(params.normalpsms) : false

activations = [hcd:'High-energy collision-induced dissociation', cid:'Collision-induced dissociation', etd:'Electron transfer dissociation']
activationtype = activations[params.activation]
massshifts = [tmt:0.0013, itraq:0.00125, false:0]
plextype = params.isobaric ? params.isobaric.replaceFirst(/[0-9]+plex/, "") : false
massshift = massshifts[plextype]
msgfprotocol = params.isobaric ? [tmt:4, itraq:2][plextype] : 0


/* PIPELINE START */

// Either feed an mzmldef file (tab separated lines with filepath\tsetname), or /path/to/\*.mzML
if (!params.mzmldef) {
Channel
  .fromPath(params.mzmls)
  .map { it -> [it, 'NA'] }
  .set { mzml_in }
} else {
Channel
  .from(file("${params.mzmldef}").readLines())
  .map { it -> it.tokenize('\t') }
  .set { mzml_in }
}

mzml_in
  .tap { sets; mzmlcounter }
  .map { it -> [file(it[0]), it[1], it[2] ? it[2] : 'NA', it[3] ? it[3].toInteger() : 'NA' ]} // create file, set plate and fraction to NA if there is none
  .tap { strips }
  .map { it -> [it[1], it[0].baseName.replaceFirst(/.*\/(\S+)\.mzML/, "\$1"), it[0], it[2], it[3]] }
  .into { mzmlfiles; groupset_mzmls; mzml_isobaric; mzml_premsgf }

mzmlcounter
  .count()
  .subscribe { println "$it mzML files in analysis" }
  .set { mzmlcount_psm }

sets
  .map{ it -> it[1] }
  .unique()
  .tap { sets_for_emtpybam; sets_for_denoms; sets_for_six }
  .collect()
  .subscribe { println "Detected setnames: ${it.join(', ')}" }


strips
  .map { it -> [it[1], it[2]] }
  .unique()
  .groupTuple()
  .set { strips_for_six }

// Get denominators if isobaric experiment
// passed in form --denoms 'set1:126:128N set2:131 set4:129N:130C:131'
if (params.isobaric && params.denoms) {
  setdenoms = [:]
  params.denoms.tokenize(' ').each{ it -> x=it.tokenize(':'); setdenoms.put(x[0], x[1..-1])}
  set_denoms = Channel.value(setdenoms)
} else if (params.isobaric) {
  setdenoms = [:]
  sets_for_denoms.reduce(setdenoms){ a, b -> a.put(b, ['_126']); return a }.set{ set_denoms }
}

if (params.pisepdb) {
  sets_for_six
    .toList()
    .map { it -> [it, normalpsms]}
    .set { normpsms }
} else {
  sets_for_six.set{ normpsms }
}

process splitSetNormalSearchPsms {
  //  normal search psm table, split on set col, collect files

  when: params.pisepdb
  input:
  set val(setnames), file('normalpsms') from normpsms
  output:
  set val(setnames), file({setnames.collect() { it + '.tsv' }}) into setnormpsms
  """
  msspsmtable split -i normalpsms --bioset
  """
}

setnormpsms
  .transpose()
  .join(strips_for_six)
  .set { setnormpsmtable } 

process splitPlateNormalSearchPsms {
  // create pep tables, split on plate, collect

  when: params.pisepdb
  input:
  set val(setname), file(normpsm), val(stripnames) from setnormpsmtable
  output:
  set val(setname), val(stripnames), file({stripnames.collect() { it + '.tsv' }}) into setplatepsms
  """
  msspsmtable split -i $normpsm --splitcol `python -c 'with open("$normpsm") as fp: h=next(fp).strip().split("\\t");print(h.index("Strip")+1)'`
  """
}

setplatepsms
  .transpose()
  .set { setplatepsmtable }

process normalSearchPsmsToPeptides {
  // create pep tables, split on plate, collect

  when: params.pisepdb
  input:
  set val(setname), val(strip), file(normpsm) from setplatepsmtable
  output:
  set val(setname), val(strip), file('peptides') into setplatepeptides
  """
  msspeptable psm2pep -i $normpsm -o peptides --scorecolpattern area --spectracol 1 
  """
}

pipep = params.pisepdb ? Channel.fromPath(params.pisepdb) : Channel.empty()
setplatepeptides
  .combine(pipep)
  .set { sixftcreation_in }

process create6FTDB {
  // create 6FT DB per peptable-plate, collect fr 

  when: params.pisepdb

  input:
  set val(setname), val(stripname), file(peptides), file(pipeptides) from sixftcreation_in

  output:
  set val(setname), val(stripname), file('target_fr*.fasta') into t_splitdb

  script:
  strip = params.strips[stripname]
  """
  pi_database_splitter.py -i $pipeptides -p $peptides --intercept $strip.intercept --width $strip.fr_width --tolerance $strip.tolerance --amount $strip.fr_amount  ${strip.reverse ? '--reverse' : ''} --deltacolpattern Delta --fdrcolpattern '^q-value' --picutoff 0.2 --fdrcutoff 0.0 --maxlen 50 --minlen 8
  """
}

// channel match plate/fr/mzML
if (params.pisepdb) {
  t_splitdb
    .transpose()
    .map { it -> ["${it[0]}_${it[1]}_${it[2].baseName.replaceFirst(/.*_fr[0]*/, "")}", it[2]]}
    .set { db_w_id }
  mzml_premsgf
    .map { it -> ["${it[0]}_${it[3]}_${it[4]}", it[0], it[1], it[2]] }  // add set_strip_fr identifier
    .into { mzml_dbid; mzml_dbfilter }
} else {
  Channel.from([['NA', tdb]]).set { db_w_id }
  mzml_premsgf
    .map { it -> ["NA", it[0], it[1], it[2]] }
    .into { mzml_dbid; mzml_dbfilter }
}

// Join DB with mzmls, use join filters out DBs without a match (in case of missing mzML fractions) 
// or duplicates (when having reruns, or when not running pI separared DBs in which case you only need 
// this process once). So we dont generate more decoys than necessary
db_w_id
  .join(mzml_dbfilter)
  .map { it -> [it[0], it[1]] }
  .set { db_filtered }

process concatFasta {
 
  input:
  set val(dbid), file(db) from db_filtered
  file knownproteins

  output:
  set val(dbid), file('td_concat.fa') into db_concatdecoy

  script:
  """
  cat $db $knownproteins > td_concat.fa
  reverse_decoy.py td_concat.fa
  cat decoy_td_concat.fa >> td_concat.fa
  rm decoy_td_concat.fa
  """
}

// Now re-match the DB (now with decoy) with mzML, this is needed to fan out if more mzMLs than
// DBs have been used so we use the cross operator
db_concatdecoy
  .cross(mzml_dbid) // gives two Arrays so unfold them in next map step
  .map { it -> [it[0][0], it[0][1], it[1][1], it[1][2], it[1][3]] } // dbid, db, set, sample, file
  .set { mzml_msgf }

process makeProtSeq {

  input:
  file knownproteins

  output:
  file('mslookup_db.sqlite') into protseqdb

  """
  msslookup protspace -i $knownproteins --minlen 8
  """
}

process makeTrypSeq {

  input:
  file knownproteins

  output:
  file('mslookup_db.sqlite') into trypseqdb

  """
  msslookup seqspace -i $knownproteins --insourcefrag
  """
}


process IsobaricQuant {

  when: !params.quantlookup && params.isobaric

  input:
  set val(setname), val(sample), file(infile), val(strip), val(fraction) from mzml_isobaric

  output:
  set val(sample), file("${infile}.consensusXML") into isobaricxml

  """
  source activate openms-blat-0.1
  IsobaricAnalyzer  -type $params.isobaric -in $infile -out "${infile}.consensusXML" -extraction:select_activation "$activationtype" -extraction:reporter_mass_shift $massshift -extraction:min_precursor_intensity 1.0 -extraction:keep_unannotated_precursor true -quantification:isotope_correction true 
  """
}

isobaricxml
  .ifEmpty(['NA', 'NA'])
  .toList()
  .flatMap { it.sort({a, b -> a[0] <=> b[0]}) }
  .map { it -> it[1] }
  .collect()
  .set { sorted_isoxml }


mzmlfiles
  .toList()
  .map { it.sort( {a, b -> a[1] <=> b[1]}) }
  .map { it -> [it.collect() { it[0] }, it.collect() { it[2] }] }
  .set{ mzmlfiles_all }


process createSpectraLookup {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: {it == 'mslookup_db.sqlite' ? 'quant_lookup.sql' : null }

  when: !params.quantlookup

  input:
  file(isobxmls) from sorted_isoxml 
  set val(setnames), file(mzmlfiles) from mzmlfiles_all
  
  output:
  file('mslookup_db.sqlite') into newspeclookup

  script:
  if(params.isobaric)
  """
  msslookup spectra -i ${mzmlfiles.join(' ')} --setnames ${setnames.join(' ')}
  msslookup isoquant --dbfile mslookup_db.sqlite -i ${isobxmls.join(' ')} --spectra ${mzmlfiles.join(' ')}
  """
  else
  """
  msslookup spectra -i ${mzmlfiles.join(' ')} --setnames ${setnames.join(' ')}
  """
}


if (!params.quantlookup) {
  newspeclookup
    .set { spec_lookup }
} else {
  Channel
    .fromPath(params.quantlookup)
    .set { spec_lookup }
} 


process msgfPlus {

  // Some versions have problems when converting to TSV, possible too long identifiers 
  // If problems arise, try to use an older version: msgf_plus:2016.10.26--py27_1

  input:
  set val(setfr_id), file(db), val(setname), val(sample), file(x) from mzml_msgf
  file mods

  output:
  set val(setname), val(sample), file("${sample}.mzid") into mzids
  set val(setname), file("${sample}.mzid"), file('out.mzid.tsv') into mzidtsvs
  
  script:
  mem = db.size() * 16 // Xmx in Bytes needed for java
  """
  msgf_plus -Xmx${task.memory.toMega()}M -d $db -s $x -o "${sample}.mzid" -thread ${task.cpus * params.threadspercore} -mod $mods -tda 0 -t 10.0ppm -ti -1,2 -m 0 -inst 3 -e 1 -protocol ${msgfprotocol} -ntt 2 -minLength 7 -maxLength 50 -minCharge 2 -maxCharge 6 -n 1 -addFeatures 1
  msgf_plus -Xmx3500M edu.ucsd.msjava.ui.MzIDToTsv -i "${sample}.mzid" -o out.mzid.tsv
  rm td_concat.c*
  """
}

mzids
  .groupTuple()
  .set { mzids_2pin }


process percolator {

  input:
  set val(setname), val(samples), file('mzid?') from mzids_2pin

  output:
  set val(setname), file('perco.xml') into percolated

  """
  echo $samples
  mkdir mzids
  count=1;for sam in ${samples.join(' ')}; do ln -s `pwd`/mzid\$count mzids/\${sam}.mzid; echo mzids/\${sam}.mzid >> metafile; ((count++));done
  msgf2pin -o percoin.xml -e trypsin -P "decoy_" metafile
  percolator -j percoin.xml -X perco.xml -N 500000 --decoy-xml-output -y
  """
}


percolated
  .tap { var_percolated }
  .set { nov_percolated }


process getVariantPercolator {

  when: params.varheaders

  input:
  set val(setname), file(x) from var_percolated

  output:
  set val(setname), val('variant'), file("${x}_h0.xml") into var_perco
  """
  msspercolator splitprotein -i $x --protheaders \'${params.varheaders}\'
  """
}


process getNovelPercolator {

  when: params.novheaders

  input:
  set val(setname), file(x) from nov_percolated

  output:
  set val(setname), val('novel'), file("${x}_h0.xml") into nov_perco
  """
  msspercolator splitprotein -i $x --protheaders \'${params.novheaders}\'
  """
}


nov_perco
  .concat(var_perco)
  .set { splitperco }


process filterPercolator {

  input:
  set val(setname), val(peptype), file(perco) from splitperco
  file 'trypseqdb' from trypseqdb
  file 'protseqdb' from protseqdb
  file knownproteins

  output:
  set val(setname), val(peptype), file('filtprot') into filtered_perco
  
  script:
  if (params.noclassfdr)
  """
  mv $perco filtprot
  """
  else
  """
  msspercolator filterseq -i $perco -o filtseq --dbfile trypseqdb --insourcefrag 2 --deamidate 
  msspercolator filterprot -i filtseq -o filtprot --fasta $knownproteins --dbfile protseqdb --minlen 8 --deamidate --enforce-tryptic
  """
}

nov_filtered_perco = Channel.create()
var_filtered_perco = Channel.create()
filtered_perco
  .choice( var_filtered_perco, nov_filtered_perco) { it -> it[1] == 'variant' ? 0 : 1 }

mzidtsvs
  .groupTuple()
  .tap { variantmzidtsv }
  .join(nov_filtered_perco)
  .set { nov_mzperco }

variantmzidtsv
  .join(var_filtered_perco)
  .concat(nov_mzperco)
  .set { allmzperco }

process svmToTSV {

  input:
  set val(setname), file('mzident?'), file('mzidtsv?'), val(peptype), file(perco) from allmzperco 

  output:
  set val(setname), val(peptype), file('mzidperco') into mzidtsv_perco

  script:
  """
#!/usr/bin/env python
from glob import glob
from lxml import etree
mzidtsvfns = sorted(glob('mzidtsv*'))
mzidfns = sorted(glob('mzident*'))
from app.readers import pycolator, xml, tsv, mzidplus
import os
ns = xml.get_namespace_from_top('$perco', None)
psms = {}
psms = {p.attrib['{%s}psm_id' % ns['xmlns']]: (float(p.find('{%s}svm_score' % ns['xmlns']).text),
            p.attrib['{%s}decoy' % ns['xmlns']] == 'true') for p in pycolator.generate_psms('$perco', ns)}
decoys = {True: 0, False: 0}
for psm in sorted([(pid, score, decoy) for pid, (score, decoy) in psms.items()], reverse=True, key=lambda x:x[1]):
    decoys[psm[2]] += 1
    try:
        psms[psm[0]] = {'decoy': psm[2], 'svm': psm[1], 'qval': decoys[True]/decoys[False]}  # T-TDC
    except ZeroDivisionError:
        psms[psm[0]] = {'decoy': psm[2], 'svm': psm[1], 'qval': 1.0}  # T-TDC
decoys = {'true': 0, 'false': 0}
for svm, decoy, psm_ids in sorted([(float(x.find('{%s}svm_score' % ns['xmlns']).text), x.attrib['{%s}decoy' % ns['xmlns']], x.find('{%s}psm_ids' % ns['xmlns'])) for x in pycolator.generate_peptides('$perco', ns)], reverse=True, key=lambda x:x[0]):
    decoys[decoy] += 1
    for pid in psm_ids:
        try:
            psms[pid.text].update({'pepqval': decoys['true']/decoys['false']})
        except ZeroDivisionError:
            psms[pid.text].update({'pepqval': 1.0})
  #      except KeyError:
  #          pass # This happens when MSGF incorrectly matches to no-ENSPs https://github.com/MSGFPlus/msgfplus/issues/78
        
oldheader = tsv.get_tsv_header(mzidtsvfns[0])
header = oldheader + ['percolator svm-score', 'PSM q-value', 'peptide q-value']
with open('mzidperco', 'w') as fp:
    fp.write('\\t'.join(header))
    for fnix, mzidfn in enumerate(mzidfns):
        mzns = mzidplus.get_mzid_namespace(mzidfn)
        siis = (sii for sir in mzidplus.mzid_spec_result_generator(mzidfn, mzns) for sii in sir.findall('{%s}SpectrumIdentificationItem' % mzns['xmlns']))
        for specidi, psm in zip(siis, tsv.generate_tsv_psms(mzidtsvfns[fnix], oldheader)):
            # percolator psm ID is: samplename_SII_scannr_rank_scannr_charge_rank
            scan, rank = specidi.attrib['id'].replace('SII_', '').split('_')
            outpsm = {k: v for k,v in psm.items()}
            spfile = os.path.splitext(psm['#SpecFile'])[0]
            try:
                percopsm = psms['{fn}_SII_{sc}_{rk}_{sc}_{ch}_{rk}'.format(fn=spfile, sc=scan, rk=rank, ch=psm['Charge'])]
            except KeyError:
                continue
            if percopsm['decoy']:
                continue
            fp.write('\\n')
            outpsm.update({'percolator svm-score': percopsm['svm'], 'PSM q-value': percopsm['qval'], 'peptide q-value': percopsm['pepqval']})
            fp.write('\\t'.join([str(outpsm[k]) for k in header]))
  """
}

mzidtsv_perco
  .combine(spec_lookup)
  .set { prepsm }

process createPSMTables {

  input:
  set val(setname), val(peptype), file('psms'), file('lookup') from prepsm
  val(mzmlcount) from mzmlcount_psm

  output:
  set val(setname), val(peptype), file("${setname}_${peptype}_psmtable.txt") into psmtable

  script:
  if(params.isobaric)
  """
  msspsmtable conffilt -i psms -o filtpsm --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'PSM q-value'
  msspsmtable conffilt -i filtpsm -o filtpep --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'peptide q-value'
  cp lookup psmlookup
  msslookup psms -i filtpep --dbfile psmlookup
  msspsmtable specdata -i filtpep --dbfile psmlookup -o prepsms.txt
  msspsmtable quant -i prepsms.txt -o ${setname}_${peptype}_psmtable.txt --dbfile psmlookup --isobaric
  sed 's/\\#SpecFile/SpectraFile/' -i ${setname}_${peptype}_psmtable.txt
  """
  else
  """
  msspsmtable conffilt -i psms -o filtpsm --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'PSM q-value'
  msspsmtable conffilt -i filtpsm -o filtpep --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'peptide q-value'
  cp lookup psmlookup
  msslookup psms -i filtpep --dbfile psmlookup
  msspsmtable specdata -i filtpep --dbfile psmlookup -o ${setname}_${peptype}_psmtable.txt
  sed 's/\\#SpecFile/SpectraFile/' -i ${setname}_${peptype}_psmtable.txt
  """
}


variantpsms = Channel.create()
novelpsms = Channel.create()
psmtable
  .tap { setmergepsmtables; peppsms }
  .choice( variantpsms, novelpsms ) { it -> it[1] == 'variant' ? 0 : 1 }

 
setmergepsmtables
  .groupTuple(by: 1)
  .set { psmmerge_in }


process mergeSetPSMtable {
  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  set val(setnames), val(peptype), file(psmtables) from psmmerge_in

  output:
  file "${peptype}_psmtable.txt" into produced_psmtables

  """
  head -n1 ${psmtables[0]} > ${peptype}_psmtable.txt
  for fn in ${psmtables.join(' ')}; do tail -n+2 \$fn >> ${peptype}_psmtable.txt; done
  """
}


process prePeptideTable {

  input:
  set val(setname), val(peptype), file('psms.txt') from peppsms 

  output:
  set val(setname), val(peptype), file('peptidetable.txt') into peptable

  script:
  """
  msspeptable psm2pep -i psms.txt -o preisoquant --scorecolpattern svm --spectracol 1 ${params.isobaric ? '--isobquantcolpattern plex' : ''}
  awk -F '\\t' 'BEGIN {OFS = FS} {print \$13,\$14,\$3,\$8,\$9,\$10,\$12,\$15,\$16,\$17,\$18,\$19,\$20,\$21,\$22,\$23}' preisoquant > preordered
  ${params.isobaric ? "msspsmtable isoratio -i psms.txt -o peptidetable.txt --targettable preordered --isobquantcolpattern plex --minint 0.1 --denompatterns ${set_denoms.value[setname].join(' ')} --protcol 12" : 'mv preordered peptidetable.txt'}
  """
}

novelpsms
  .into{novelpsmsFastaBedGFF; novelpsms_specai}

novelprepep = Channel.create()
presai_peptable = Channel.create()
peptable
  .choice( presai_peptable, novelprepep ) { it -> it[1] == 'variant' ? 0 : 1 }
novelprepep
  .join(novelpsmsFastaBedGFF)
  .set { novelFaBdGfPep }

process createFastaBedGFF {
 
  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { it == "${setname}_novel_peptides.gff3" ? "${setname}_novel_peptides.gff3" : null}
 
  input:
  set val(setname), val(peptype), file(peptides) , val(psmtype), file(psms) from novelFaBdGfPep
  file gtffile
  file tdb from tdb
 
  output:
  set val(setname), file('novel_peptides.fa') into novelfasta
  set val(setname), file('novel_peptides.bed') into novelbed
  set val(setname), file("${setname}_novel_peptides.gff3") into novelGFF3
  set val(setname), file('novel_peptides.tab.txt') into novelpep
  set val(setname), file('novpep_perco_quant.txt') into novelpep_percoquant
 
  """
  map_novelpeptide2genome.py --input $psms --gtf $gtffile --fastadb $tdb --tab_out novel_peptides.tab.txt --fasta_out novel_peptides.fa --gff3_out ${setname}_novel_peptides.gff3 --bed_out novel_peptides.bed
  sort -k 1b,1 <(tail -n+2 $peptides) |cut -f 1,14-500 > peptable_sorted
  sort -k 2b,2 <(tail -n+2 novel_peptides.tab.txt) > novpep_sorted
  paste <(cut -f 2 novpep_sorted) <(cut -f1,3-500 novpep_sorted) > novpep_pepcols
  join novpep_pepcols peptable_sorted -j 1 -a1 -o auto -e 'NA' -t \$'\\t' > novpep_pqjoin
  # Cut only bare peptide col and q-values/isoquant
  paste <(cut -f 2 novpep_pqjoin) <(cut -f8-500 novpep_pqjoin) > novpep_joined_pepcols
  paste <(head -n1 novel_peptides.tab.txt | cut -f1)  <(cut -f 14-500 $peptides |head -n1) > header
  cat header novpep_joined_pepcols > novpep_perco_quant.txt
  """
}

novelpep
  .into {blastnovelpep; blatnovelpep; annonovelpep; snpnovelpep}
novelfasta
  .into {blastnovelfasta; blatnovelfasta}

process BlastPNovel {

  input:
  set val(setname), file(novelfasta) from blastnovelfasta
  file blastdb

  output:
  set val(setname), file('blastp_out.txt') into novelblast
  
  """
  makeblastdb -in $blastdb -dbtype prot
  blastp -db $blastdb -query $novelfasta -outfmt '6 qseqid sseqid pident qlen slen qstart qend sstart send mismatch positive gapopen gaps qseq sseq evalue bitscore' -num_threads 4 -max_target_seqs 1 -evalue 1000 -out blastp_out.txt
  """
}

novelpsms_specai
  .map { it -> [it[0], it[2]] }
  .join(blastnovelpep)
  .join(novelblast)
  .set { novelblastout }

process ParseBlastpOut {
 
 input:
 set val(setname), file(psms), file(novelpep), file(novelblast) from novelblastout
 file blastdb

 output:
 set val(setname), file('peptable_blastp.txt') into peptable_blastp
 set val(setname), file('single_mismatch_novpeps.txt') into novpeps_singlemis

 """
 parse_BLASTP_out.py --input $novelpep --blastp_result $novelblast --fasta $blastdb --output peptable_blastp.txt
 extract_1mismatch_novpsm.py peptable_blastp.txt $psms single_mismatch_novpeps.txt
 """

}

groupset_mzmls
  .map { it -> [it[0], it[1], it[2]] } // strip fraction, strip
  .groupTuple()
  .tap { var_specaimzmls }
  .join(novpeps_singlemis)
  .set { grouped_saavnov_mzml_peps }
  

process ValidateSingleMismatchNovpeps {
  
  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { it == "precursorError.histogram.plot.pdf" ? "${setname}_novel_precursorError_plot.pdf" : it }
  
  input:
  set val(setname), val(samples), file(mzmls), file(peps) from grouped_saavnov_mzml_peps

  output:
  set val(setname), file("${setname}_novel_saav_specai.txt") into singlemis_specai
  file 'precursorError.histogram.plot.pdf' optional true into novel_specai_plot

  """
  mkdir mzmls
  for fn in $mzmls; do ln -s `pwd`/\$fn mzmls/; done
  Rscript /SpectrumAI/SpectrumAI.R mzmls $peps ${setname}_novel_saav_specai.txt || cp $peps singlemis_specai.txt
  """
}

singlemis_specai
  .join(peptable_blastp)
  .set { nov_specaiparse }

process novpepSpecAIOutParse {

  input:
  set val(setname), file(x), file('peptide_table.txt') from nov_specaiparse
  
  output:
  set val(setname), file('novpep_specai.txt') into novpep_singlemisspecai

  """
  parse_spectrumAI_out.py --spectrumAI_out $x --input peptide_table.txt --output novpep_sa
  cut -f 1,8-19 novpep_sa > novpep_specai.txt
  """
}


process BLATNovel {

  input:
  set val(setname), file(novelfasta) from blatnovelfasta
  file genomefa

  output:
  set val(setname), file('blat_out.pslx') into novelblat

  """
  source activate openms-blat-0.1
  blat $genomefa $novelfasta -t=dnax -q=prot -tileSize=5 -minIdentity=99 -out=pslx blat_out.pslx 
  """
}

novelblat
  .join(blatnovelpep)
  .set { novblatparse }

process parseBLATout {

 input:
 set val(setname), file(novelblat), file(novelpep) from novblatparse

 output:
 set val(setname), file('peptable_blat.txt') into peptable_blat

 """
 parse_BLAT_out.py $novelblat $novelpep peptable_blat.txt

 """
}

process labelnsSNP {
  
  input:
  set val(setname), file(peptable) from snpnovelpep
  file snpfa

  output:
  set val(setname), file('nssnp.txt') into ns_snp_out

  """
  label_nsSNP_pep.py --input $peptable --nsSNPdb $snpfa --output nssnp.txt
  """
}

bwfile = Channel.fromPath(params.bigwigs)
novelGFF3
  .combine(bwfile)
  .into { novelGFF3_phast; novelGFF3_phylo; novelGFF3_bams }

process phastcons {
  
  input:
  set val(setname), file(novelgff), file('bigwigs') from novelGFF3_phast

  output:
  set val(setname), file ('phastcons.txt') into phastcons_out

  """
  calculate_phastcons.py $novelgff bigwigs/hg19.100way.phastCons.bw phastcons.txt
  """
}

process phyloCSF {
  
  input:
  set val(setname), file(novelgff), file('bigwigs') from novelGFF3_phylo

  output:
  set val(setname), file('phylocsf.txt') into phylocsf_out

  """
  calculate_phylocsf.py $novelgff bigwigs phylocsf.txt
  """

}


bamFiles = params.bamfiles ? Channel.fromPath(params.bamfiles).map { fn -> [ fn, fn + '.bai' ] } : Channel.empty()

process scanBams {

  when: params.bamfiles

  input:
  set val(setname), file(gff) from novelGFF3_bams
  file bams from bamFiles.collect()
  
  output:
  set val(setname), file('scannedbams.txt') into scannedbams

  """
  ls *.bam > bamfiles.txt
  scan_bams.py  --gff_input $gff --bam_files bamfiles.txt --output scannedbams.txt
  """
}

annoperl = Channel.fromPath("$params.annovar_dir/annotate_variation.pl")
annohumdb = Channel.fromPath("$params.annovar_dir/humandb/")

novelbed
  .combine(annoperl)
  .combine(annohumdb)
  .set { anno_in }

process annovar {
  
  input:
  set val(setname), file(novelbed), file(perlscript), file(humdb) from anno_in

  output:
  set val(setname), file('novpep_annovar.variant_function') into annovar_out

  """
  ./annotate_variation.pl -out novpep_annovar -build hg19 $novelbed humandb/
  """

}

annovar_out
  .join(annonovelpep)
  .set { parseanno }

process parseAnnovarOut {
  
  input:
  set val(setname), file(anno), file(novelpep) from parseanno

  output:
  set val(setname), file('parsed_annovar.txt') into annovar_parsed

  """
  parse_annovar_out.py --input $novelpep --output parsed_annovar.txt --annovar_out $anno 
  """
}


ns_snp_out
  .join(novpep_singlemisspecai)
  .join(peptable_blat)
  .join(annovar_parsed)
  .join(phastcons_out)
  .join(phylocsf_out)
  .join(novelpep_percoquant)
  .set { combined_novelprebam }

combined_novel = (params.bamfiles ? combined_novelprebam.join(scannedbams) : combined_novelprebam)


process combineResults{
  
  publishDir "${params.outdir}", mode: 'copy', overwrite: true
  
  input:
  set val(setname), file(a), file(b), file(c), file(d), file(e), file(f), file(g), file(h) from combined_novel
  
  output:
  set val('nov'), val(setname), file("${setname}_novel_peptides.txt") into novpeps_finished 
  
  script:
  if (!params.bamfiles)
  """
  for fn in $a $b $c $d $e $f $g; do sort -k 1b,1 \$fn > tmpfn; mv tmpfn \$fn; done
  join $a $b -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined1
  join joined1 $c -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined2
  join joined2 $d -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined3
  join joined3 $e -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined4
  join joined4 $f -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined5
  join joined5 $g -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined6
  grep '^Bare peptide' joined6 > ${setname}_novel_peptides.txt
  grep -v '^Bare peptide' joined6 >> ${setname}_novel_peptides.txt
  """

  else
  """
  for fn in $a $b $c $d $e $f $g $h; do sort -k 1b,1 \$fn > tmpfn; mv tmpfn \$fn; done
  join $a $b -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined1
  join joined1 $c -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined2
  join joined2 $d -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined3
  join joined3 $e -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined4
  join joined4 $f -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined5
  join joined5 $g -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined6
  join joined6 $h -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined7
  grep '^Bare peptide' joined7 > ${setname}_novel_peptides.txt
  grep -v '^Bare peptide' joined7 >> ${setname}_novel_peptides.txt
  """
}

process prepSpectrumAI {

  input:
  set val(setname), val(peptype), file(psms) from variantpsms
  
  output:
  set val(setname), file('specai_in.txt') into var_specai_input
  
  script:
  if (params.saavheader)
  """
  cat <(head -n1 $psms) <(grep $params.saavheader $psms) > saavpsms
  label_sub_pos.py --input_psm saavpsms --output specai_in.txt ${params.splitchar ? "--splitchar ${params.splitchar}" : ''}
  """
  else
  """
  label_sub_pos.py --input_psm $psms --output specai_in.txt
  """
}


var_specaimzmls
  .join(var_specai_input)
  .set { var_specai_inmzml }

process SpectrumAI {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { it == "precursorError.histogram.plot.pdf" ? "${setname}_variant_precursorError_plot.pdf" : it }

  input:
  set val(setname), val(samples), file(mzmls), file(specai_in) from var_specai_inmzml

  output:
  set val(setname), file("${setname}_variant_specairesult.txt") into specai
  file "precursorError.histogram.plot.pdf" into specai_plot

  """
  mkdir mzmls
  for fn in $mzmls; do ln -s `pwd`/\$fn mzmls/; done
  ls mzmls
  Rscript /SpectrumAI/SpectrumAI.R mzmls $specai_in ${setname}_variant_specairesult.txt
  """
}

specai
  .join(presai_peptable)
  .set { specai_peptable }

process mapVariantPeptidesToGenome {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  set val(setname), file(x), val(peptype), file(peptides) from specai_peptable
  file cosmic
  file dbsnp
  
  output:
  set val('var'), val(setname), file("${setname}_variant_peptides.txt") into varpeps_finished
  file "${setname}_variant_peptides.saav.pep.hg19cor.vcf" into saavvcfs_finished

  """
  ${params.saavheader ? "cat <(head -n1 ${peptides}) <(grep ${params.saavheader} ${peptides}) > saavpeps" : "mv ${peptides} saavpeps" }
  parse_spectrumAI_out.py --spectrumAI_out $x --input saavpeps --output setsaavs
  ${params.saavheader ? "cat setsaavs <(grep -v ${params.saavheader} ${peptides} | sed \$'s/\$/\tNA/') > ${setname}_variant_peptides.txt" : "mv setsaavs ${setname}_variant_peptides.txt"}
  map_cosmic_snp_tohg19.py --input ${setname}_variant_peptides.txt --output ${setname}_variant_peptides.saav.pep.hg19cor.vcf --cosmic_input $cosmic --dbsnp_input $dbsnp
  # Remove PSM-table specific stuff (RT, precursor, etc etc) from variant PEPTIDE table
  cut -f 1,2,14-5000 ${setname}_variant_peptides.txt > pepsfix
  mv pepsfix ${setname}_variant_peptides.txt
  """
}

novpeps_finished
  .concat(varpeps_finished) 
  .groupTuple()
  .set { setmerge_peps }


accession_keymap = ['var': 'Peptide sequence', 'nov': 'Peptide']
acc_removemap = ['nov': 'Bare peptide', 'var': 'Mod.peptide']


process mergeSetPeptidetable {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  set val(peptype), val(setnames), file('peps?') from setmerge_peps

  output:
  file "${peptype}_peptidetable.txt" into produced_peptables

  """
  # build non-changing fields (seq based fields) table:
  fixfields=`head -n1 peps1 |tr -s '\\t' '\\n' | egrep -vn '(Setname|Spectrum|q-val|plex|${acc_removemap[peptype]})' | cut -f 1 -d ':'`
  fixfields=`echo \$fixfields | sed 's/ /,/g'`
  head -n1 peps1 | cut -f `echo \$fixfields` > fixheader
  count=1; for setn in ${setnames.join(' ')} ; do
  cut -f  `echo \$fixfields` peps\$count | tail -n+2 >> fixpeps
  ((count++))
  done
  if [ ${peptype} == 'nov' ]
  then
     cat fixheader <(sort -u -k1b,1 fixpeps) > temp
     group_novpepToLoci.py  --input temp --output temp.loci --distance 10kb
     head -n1 temp.loci > fixheader
     tail -n+2 temp.loci > fixpeps
  fi
  sort -u -k1b,1 fixpeps > temp
  mv temp fixpeps

  ## Build changing fields table
  touch peptable
  count=1; for setn in ${setnames.join(' ')}; do
    varfields=`head -n1 peps\$count |tr -s '\\t' '\\n' | egrep -n '(${accession_keymap[peptype]}|Spectrum|q-val|plex)' | cut -f 1 -d ':'`
    varfields=`echo \$varfields| sed 's/ /,/g'`
    # first add to header, cut from f2 to remove join-key pep seq field
    head -n1 peps\$count | cut -f `echo \$varfields` | cut -f 2-5000| sed "s/^\\(\\w\\)/\${setn}_\\1/;s/\\(\\s\\)/\\1\${setn}_/g" > varhead
    paste fixheader varhead > newheader && mv newheader fixheader
    # then join the values
    tail -n+2 peps\$count | cut -f `echo \$varfields` | sort -k1b,1 > sortpep; join peptable sortpep -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined
    mv joined peptable
    ((count++))
  done
  join fixpeps peptable -a1 -a2 -o auto -e 'NA' -t \$'\\t' > fixvarpeps
  cat fixheader fixvarpeps > ${peptype}_peptidetable.txt
  """
}


produced_psmtables
  .concat(produced_peptables)
  .subscribe { println "Pipeline output ready: ${it}" }
