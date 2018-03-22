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

FIXME:
 - make sure docker files are ok and in right place etc. Try them!
*/

nf_required_version = '0.26.0'
if( ! nextflow.version.matches(">= ${nf_required_version}") ){
  println("Nextflow version too old, ${nf_required_version} required")
  exit(1)
}


/* SET DEFAULT PARAMS */
mods = file('Mods.txt')
params.ppoolsize = 8
params.isobaric = false
params.activation = 'hcd'
params.bamfiles = false

knownproteins = file(params.knownproteins)
blastdb = file(params.blastdb)
gtffile = file(params.gtf)
snpfa = file(params.snpfa)
dbsnp = file(params.dbsnp)
cosmic = file(params.cosmic)
genomefa = file(params.genome)
tdb = file(params.tdb)

activations = [hcd:'High-energy collision-induced dissociation', cid:'Collision-induced dissociation', etd:'Electron transfer dissociation']
activationtype = activations[params.activation]
massshifts = [tmt:0.0013, itraq:0.00125, false:0]
plextype = params.isobaric ? params.isobaric.replaceFirst(/[0-9]+plex/, "") : false
massshift = massshifts[plextype]
msgfprotocol = [tmt:4, itraq:2, false:0][plextype]

/* PIPELINE START */
Channel
  .fromPath(params.mzmls)
  .count()
  .set{ amount_mzml }


process concatFasta {
 
  container 'ubuntu:latest'

  input:
  file tdb
  file knownproteins

  output:
  file('db.fa') into targetdb

  script:
  """
  cat $tdb $knownproteins > db.fa
  """
}


process makeDecoyReverseDB {
  container 'biopython/biopython'

  input:
  file db from targetdb

  output:
  file('concatdb.fasta') into concatdb

  """
  #!/usr/bin/env python3
  from Bio import SeqIO
  with open('$db') as fp, open('concatdb.fasta', 'w') as wfp:
    for target in SeqIO.parse(fp, 'fasta'):
      SeqIO.write(target, wfp, 'fasta')
      decoy = target[::-1] 
      decoy.description = decoy.description.replace('ENS', 'decoy_ENS')
      decoy.id = 'decoy_{}'.format(decoy.id)
      SeqIO.write(decoy, wfp, 'fasta')
  """
}


Channel
  .fromPath(params.mzmls)
  .map { it -> [it.baseName.replaceFirst(/.*fr(\d\d).*/, "\$1").toInteger(), it.baseName.replaceFirst(/.*\/(\S+)\.mzML/, "\$1"), it] }
  .into{ mzmlfiles; mzml_isobaric; mzml_msgf }

mzmlfiles
  .buffer(size: amount_mzml.value)
  .flatMap { it.sort( {a, b -> a[1] <=> b[1]}) }
  .map { it -> it[2] }
  .collect()
  .into { mzmlfiles_all; specaimzmls; singlemismatch_nov_mzmls }


process makeProtSeq {

  container 'quay.io/biocontainers/msstitch:2.5--py36_0'

  input:
  file knownproteins

  output:
  file('mslookup_db.sqlite') into protseqdb

  """
  msslookup protspace -i $knownproteins --minlen 8
  """
}

process makeTrypSeq {

  container 'quay.io/biocontainers/msstitch:2.5--py36_0'

  input:
  file knownproteins

  output:
  file('mslookup_db.sqlite') into trypseqdb

  """
  msslookup seqspace -i $knownproteins --insourcefrag
  """
}


process IsobaricQuant {

  container 'quay.io/biocontainers/openms:2.2.0--py27_boost1.64_0'

  when: params.isobaric

  input:
  set val(fr), val(sample), file(infile) from mzml_isobaric

  output:
  set val(sample), file("${infile}.consensusXML") into isobaricxml

  """
  IsobaricAnalyzer  -type $params.isobaric -in $infile -out "${infile}.consensusXML" -extraction:select_activation "$activationtype" -extraction:reporter_mass_shift $massshift -extraction:min_precursor_intensity 1.0 -extraction:keep_unannotated_precursor true -quantification:isotope_correction true 
  """
}

isobaricamount = params.isobaric ? amount_mzml.value : 1

isobaricxml
  .ifEmpty(['NA', 'NA'])
  .buffer(size: isobaricamount)
  .flatMap { it.sort({a, b -> a[0] <=> b[0]}) }
  .map { it -> it[1] }
  .collect()
  .set { sorted_isoxml }


process createSpectraLookup {

  container 'quay.io/biocontainers/msstitch:2.5--py36_0'

  input:
  file(isobxmls) from sorted_isoxml 
  file(mzmlfiles) from mzmlfiles_all
  
  output:
  file 'mslookup_db.sqlite' into spec_lookup

  script:
  if(params.isobaric)
  """
  msslookup spectra -i ${mzmlfiles.join(' ')} --setnames  ${['setA'].multiply(amount_mzml.value).join(' ')}
  msslookup isoquant --dbfile mslookup_db.sqlite -i ${isobxmls.join(' ')} --spectra ${mzmlfiles.join(' ')}
  """
  else
  """
  msslookup spectra -i ${mzmlfiles_all.join(' ')} --setnames  ${['setA'].multiply(amount_mzml.value).join(' ')}
  """
}


process msgfPlus {

  /* Latest version has problems when converting to TSV, possible too long identifiers 
     So we use an older version
     LATEST TESTED: container 'quay.io/biocontainers/msgf_plus:2017.07.21--py27_0'
  */
  container 'quay.io/biocontainers/msgf_plus:2016.10.26--py27_1'

  input:
  set val(fraction), val(sample), file(x) from mzml_msgf
  file(db) from concatdb
  file mods

  output:
  set val(fraction), val(sample), file("${sample}.mzid") into mzids
  set val(sample), file("${sample}.mzid"), file('out.mzid.tsv') into mzidtsvs
  
  """
  msgf_plus -Xmx16G -d $db -s $x -o "${sample}.mzid" -thread 12 -mod $mods -tda 0 -t 10.0ppm -ti -1,2 -m 0 -inst 3 -e 1 -protocol ${msgfprotocol} -ntt 2 -minLength 7 -maxLength 50 -minCharge 2 -maxCharge 6 -n 1 -addFeatures 1
  msgf_plus -Xmx3500M edu.ucsd.msjava.ui.MzIDToTsv -i "${sample}.mzid" -o out.mzid.tsv
  """
}

mzids
  .map { it -> [it[1], it[2]] }
  .tap { mzids_perco }
  .buffer(size: amount_mzml.value)
  .map { it.sort( {a, b -> a[0] <=> b[0]}) }
  .map { it -> [it.collect() { it[0] }, it.collect() { it[1] }] }
  .set { mzids_2pin }

/*
tmzids = Channel.create()
dmzids = Channel.create()

tmzids
  .buffer(size: amount_mzml.value)
  .flatMap { it.sort( {a, b -> a['sample'] <=> b['sample']}) }
  .buffer(size: params.ppoolsize, remainder: true) 
  .map { it -> [it.collect() { it['fn'] }, it.collect() { it['sample'] }] }
  .set { buffer_mzid_target }
dmzids
  .buffer(size: amount_mzml.value)
  .flatMap { it.sort({a, b -> a['sample'] <=> b['sample'] }) }
  .buffer(size: params.ppoolsize, remainder: true)
  .map { it -> [it.collect() { it['fn'] }, it.collect() { it['sample'] }] }
  .set { buffer_mzid_decoy }
*/


process percolator {

  container 'quay.io/biocontainers/percolator:3.1--boost_1.623'

  input:
  set val(samples), file('mzid?') from mzids_2pin

  output:
  file('perco.xml') into percolated

  """
  echo $samples
  mkdir mzids
  count=1;for sam in ${samples.join(' ')}; do ln -s `pwd`/mzid\$count mzids/\${sam}.mzid; echo mzids/\${sam}.mzid >> metafile; ((count++));done
  msgf2pin -o percoin.xml -e trypsin -P "decoy_" metafile
  percolator -j percoin.xml -X perco.xml -N 500000 --decoy-xml-output -y
  """
}

process filterPercolator {

  container 'quay.io/biocontainers/msstitch:2.5--py36_0'

  input:
  file x from percolated
  file 'trypseqdb' from trypseqdb
  file 'protseqdb' from protseqdb
  file knownproteins

  output:
  set val('novel'), file('fp_th0.xml') into var_filtered_perco
  set val('variant'), file('fp_th1.xml') into nov_filtered_perco
  """
  msspercolator splitprotein -i perco.xml --protheaders '^PGOHUM;^lnc;^decoy_PGOHUM;^decoy_lnc' '^COSMIC;^CanProVar;^decoy_COSMIC;^decoy_CanProVar'
  msspercolator filterseq -i perco.xml_h0.xml -o fs_th0.xml --dbfile trypseqdb --insourcefrag 2 --deamidate 
  msspercolator filterseq -i perco.xml_h1.xml -o fs_th1.xml --dbfile trypseqdb --insourcefrag 2 --deamidate 
  msspercolator filterprot -i fs_th0.xml -o fp_th0.xml --fasta $knownproteins --dbfile protseqdb --minlen 8 --deamidate --enforce-tryptic
  msspercolator filterprot -i fs_th1.xml -o fp_th1.xml --fasta $knownproteins --dbfile protseqdb --minlen 8 --deamidate --enforce-tryptic
  """
}

var_filtered_perco
  .concat(nov_filtered_perco)
  .set { filtered_perco }

  //set val(sample), (file('out.mzid.tsv') into mzidtsvs
mzidtsvs
  .buffer(size: amount_mzml.value)
  .map { it -> [it.collect() { it[1] }, it.collect() { it[2] }] }
  .combine(filtered_perco)
  .set { allmzidtsv }


process svmToTSV {

  container 'quay.io/biocontainers/msstitch:2.5--py36_0'

  input:
  set file('mzident?'), file('mzidtsv?'), val(peptype), file(perco) from allmzidtsv

  output:
  set val(peptype), file('mzidperco') into mzidtsv_perco

  script:
  """
#!/usr/bin/env python
from glob import glob
mzidtsvfns = sorted(glob('mzidtsv*'))
mzidfns = sorted(glob('mzident*'))
from app.readers import pycolator, xml, tsv, mzidplus
import os
ns = xml.get_namespace_from_top('$perco', None) 
psms = {p.attrib['{%s}psm_id' % ns['xmlns']]: p for p in pycolator.generate_psms('$perco', ns)}
decoys = {True: 0, False: 0}
for psm in sorted([(pid, float(p.find('{%s}svm_score' % ns['xmlns']).text), p) for pid, p in psms.items()], reverse=True, key=lambda x:x[1]):
    pdecoy = psm[2].attrib['{%s}decoy' % ns['xmlns']] == 'true'
    decoys[pdecoy] += 1
    psms[psm[0]] = {'decoy': pdecoy, 'svm': psm[1], 'qval': decoys[True]/decoys[False]}  # T-TDC
decoys = {'true': 0, 'false': 0}
for svm, pep in sorted([(float(x.find('{%s}svm_score' % ns['xmlns']).text), x) for x in pycolator.generate_peptides('$perco', ns)], reverse=True, key=lambda x:x[0]):
    decoys[pep.attrib['{%s}decoy' % ns['xmlns']]] += 1
    [psms[pid.text].update({'pepqval': decoys['true']/decoys['false']}) for pid in pep.find('{%s}psm_ids' % ns['xmlns'])]
oldheader = tsv.get_tsv_header(mzidtsvfns[0])
header = oldheader + ['percolator svm-score', 'PSM q-value', 'peptide q-value']
with open('mzidperco', 'w') as fp:
    fp.write('\\t'.join(header))
    for fnix, mzidfn in enumerate(mzidfns):
        mzns = mzidplus.get_mzid_namespace(mzidfn)
        siis = (sii for sir in mzidplus.mzid_spec_result_generator(mzidfn, mzns) for sii in sir.findall('{%s}SpectrumIdentificationItem' % mzns['xmlns']))
        for specidi, psm in zip(siis, tsv.generate_tsv_psms(mzidtsvfns[fnix], oldheader)):
            # percolator psm ID is: samplename_SII_scannr_rank_scannr_charge_rank
            print(specidi)
            print(psm)
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

process createPSMPeptideTable {

  container 'quay.io/biocontainers/msstitch:2.5--py36_0'
  
  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  set val(peptype), file('psms') from mzidtsv_perco
  file 'lookup' from spec_lookup

  output:
  set val(peptype), file("${peptype}_psmtable.txt") into psmtable

  script:
  if(params.isobaric)
  """
  msspsmtable conffilt -i psms -o filtpsm --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'PSM q-value'
  msspsmtable conffilt -i filtpsm -o filtpep --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'peptide q-value'
  cp lookup psmlookup
  msslookup psms -i filtpep --dbfile psmlookup
  msspsmtable specdata -i filtpep --dbfile psmlookup -o prepsms.txt
  msspsmtable quant -i prepsms.txt -o ${peptype}_psmtable.txt --dbfile psmlookup --isobaric
  sed 's/\\#SpecFile/SpectraFile/' -i ${peptype}_psmtable.txt
  """
  else
  """
  msspsmtable conffilt -i psms -o filtpsm --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'PSM q-value'
  msspsmtable conffilt -i filtpsm -o filtpep --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'peptide q-value'
  cp lookup psmlookup
  msslookup psms -i filtpep --dbfile psmlookup
  msspsmtable specdata -i filtpep --dbfile psmlookup -o ${peptype}_psmtable.txt
  sed 's/\\#SpecFile/SpectraFile/' -i ${peptype}_psmtable.txt
  """
}

variantpsms = Channel.create()
novelpsms = Channel.create()
psmtable
  .tap { both_psmtables }
  .choice( variantpsms, novelpsms ) { it -> it[0] == 'variant' ? 0 : 1 }
 both_psmtables
  .map { it -> it[1] }
  .collect()
  .set { psms_prepep }

 
process prePeptideTable {

  container 'quay.io/biocontainers/msstitch:2.5--py36_0'
  
  input:
  file('psms?') from psms_prepep

  output:
  file('prepeptidetable.txt') into prepeptable

  script:
  if(params.isobaric)
  """
  msspsmtable merge -o psms.txt -i psms* 
  msspeptable psm2pep -i psms.txt -o prepeptidetable.txt --scorecolpattern svm --spectracol 1 --isobquantcolpattern plex
  """
  else
  """
  msspsmtable merge -o psms.txt -i psms* 
  msspeptable psm2pep -i psms.txt -o prepeptidetable.txt --scorecolpattern svm --spectracol 1
  """
}


process createPeptideTable{

  container 'ubuntu:latest'

  input:
  file 'prepeptidetable.txt' from prepeptable

  output:
  file 'peptide_table.txt' into peptable

  """
  paste <( cut -f 12 prepeptidetable.txt) <( cut -f 13 prepeptidetable.txt) <( cut -f 3,7-9,11,14-22 prepeptidetable.txt) > peptide_table.txt
  """
}


novelpsms
  .map { it -> it[1] }
  .into{novelpsmsFastaBedGFF; novelpsms_specai}


process createFastaBedGFF {
  container 'pgpython'
 
  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { it == "novel_peptides.gff3" ? "novel_peptides.gff3" : null}
 
  input:
  file novelpsmsFastaBedGFF
  file gtffile
  file tdb
 
  output:
  file 'novel_peptides.fa' into novelfasta
  file 'novel_peptides.bed' into novelbed
  file 'novel_peptides.gff3' into novelGFF3
  file 'novel_peptides.tab.txt' into novelpep
 
  """
  python3 /pgpython/map_novelpeptide2genome.py --input $novelpsmsFastaBedGFF --gtf $gtffile --fastadb $tdb --tab_out novel_peptides.tab.txt --fasta_out novel_peptides.fa --gff3_out novel_peptides.gff3 --bed_out novel_peptides.bed
  """
}

novelpep
  .into {blastnovelpep; blatnovelpep; annonovelpep; snpnovelpep}
novelfasta
  .into {blastnovelfasta; blatnovelfasta}

process BlastPNovel {

  container 'quay.io/biocontainers/blast:2.7.1--boost1.64_1'

  input:
  file novelfasta from blastnovelfasta
  file blastdb

  output:
  file 'blastp_out.txt' into novelblast
  
  """
  makeblastdb -in $blastdb -dbtype prot
  blastp -db $blastdb -query $novelfasta -outfmt '6 qseqid sseqid pident qlen slen qstart qend sstart send mismatch positive gapopen gaps qseq sseq evalue bitscore' -num_threads 8 -max_target_seqs 1 -evalue 1000 -out blastp_out.txt
  """
}

process ParseBlastpOut {
 container 'pgpython'
 
 input:
 file novelpsms from novelpsms_specai
 file novelpep from blastnovelpep
 file novelblast from novelblast
 file blastdb

 output:
 file 'peptable_blastp.txt' into peptable_blastp
 file 'single_mismatch_novpeps.txt' into novpeps_singlemis

 """
 python3 /pgpython/parse_BLASTP_out.py --input $novelpep --blastp_result $novelblast --fasta $blastdb --output peptable_blastp.txt
 python3 /pgpython/extract_1mismatch_novpsm.py peptable_blastp.txt $novelpsms single_mismatch_novpeps.txt
 """

}

process ValidateSingleMismatchNovpeps {
  container 'spectrumai'
  
  input:
  file x from novpeps_singlemis
  file mzml from singlemismatch_nov_mzmls

  output:
  file 'singlemis_specai.txt' into singlemis_specai

  """
  mkdir mzmls
  for fn in $mzml; do ln -s `pwd`/\$fn mzmls/; done
  Rscript /SpectrumAI/SpectrumAI.R mzmls $x singlemis_specai.txt || cp $x singlemis_specai.txt
  """
}

process novpepSpecAIOutParse {
  container 'pgpython'

  input:
  file x from singlemis_specai 
  file 'peptide_table.txt' from peptable_blastp 
  
  output:
  file 'novpep_specai.txt' into novpep_singlemisspecai

  """
  python3 /pgpython/parse_spectrumAI_out.py --spectrumAI_out $x --input peptide_table.txt --output novpep_specai.txt
  """
}

process BLATNovel {
  container 'quay.io/biocontainers/blat:35--1'

  input:
  file novelfasta from blatnovelfasta
  file genomefa

  output:
  file 'blat_out.pslx' into novelblat

  """
  blat $genomefa $novelfasta -t=dnax -q=prot -tileSize=5 -minIdentity=99 -out=pslx blat_out.pslx 
  """
}

process parseBLATout {
 container 'pgpython'

 input:
 file novelblat from novelblat
 file novelpep from blatnovelpep

 output:
 file 'peptable_blat.txt' into peptable_blat

 """
 python3 /pgpython/parse_BLAT_out.py $novelblat $novelpep peptable_blat.txt

 """
}

process labelnsSNP {
  
  container 'pgpython'
  
  input:
  file peptable from snpnovelpep
  file snpfa

  output:
  file 'nssnp.txt' into ns_snp_out

  """
  python3 /pgpython/label_nsSNP_pep.py --input $peptable --nsSNPdb $snpfa --output nssnp.txt
  """
}

novelGFF3
  .into { novelGFF3_phast; novelGFF3_phylo; novelGFF3_bams }

process phastcons {
  container 'pgpython'
  
  input:
  file novelgff from novelGFF3_phast
  output:
  file 'phastcons.txt' into phastcons_out

  """
  python3 /pgpython/calculate_phastcons.py $novelgff /bigwigs/hg19.100way.phastCons.bw phastcons.txt
  """
}

process phyloCSF {
  
  container 'pgpython'

  input:
  file novelgff from novelGFF3_phylo

  output:
  file 'phylocsf.txt' into phylocsf_out

  """
  python3 /pgpython/calculate_phylocsf.py $novelgff /bigwigs phylocsf.txt
  """

}


if (params.bamfiles) {
  bamFiles = Channel
    .fromPath(params.bamfiles)
    .map { fn -> [ fn, fn + '.bai' ] }
    .collect()
} else {
  bamFiles = Channel.empty()
}


process scanBams {
  container 'pgpython'

  when: params.bamfiles

  input:
  file gff from novelGFF3_bams
  file bams from bamFiles
  
  output:
  file 'scannedbams.txt' into scannedbams

  """
  ls *.bam > bamfiles.txt
  python3 /pgpython/scan_bams.py  --gff_input $gff --bam_files bamfiles.txt --output scannedbams.txt
  """
}


process annovar {
  
  container 'annovar'
  
  input:
  file novelbed
  output:
  file 'novpep_annovar.variant_function' into annovar_out

  """
  /annovar/annotate_variation.pl -out novpep_annovar -build hg19 $novelbed /annovar/humandb/
  """

}

process parseAnnovarOut {
  
  container 'pgpython'
  
  input:
  file anno from annovar_out
  file novelpep from annonovelpep

  output:
  file 'parsed_annovar.txt' into annovar_parsed

  """
  python3 /pgpython/parse_annovar_out.py --input $novelpep --output parsed_annovar.txt --annovar_out $anno 
  """
}

scannedbams
  .ifEmpty('No bams')
  .set { bamsOrEmpty }

process combineResults{
  
  container 'pgpython'

  input:
  file a from ns_snp_out
  file b from novpep_singlemisspecai
  file c from peptable_blat
  file d from annovar_parsed
  file e from phastcons_out
  file f from phylocsf_out
  file g from bamsOrEmpty
  
  output:
  file 'combined' into combined_novelpep_output
  
  script:
  if (!params.bamfiles)
  """
  for fn in $a $b $c $d $e $f; do sort -k 1b,1 \$fn > tmpfn; mv tmpfn \$fn; done
  join $a $b -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined1
  join joined1 $c -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined2
  join joined2 $d -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined3
  join joined3 $e -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined4
  join joined4 $f -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined5
  grep '^Peptide' joined5 > combined
  grep -v '^Peptide' joined5 >> combined
  """

  else
  """
  for fn in $a $b $c $d $e $f $g; do sort -k 1b,1 \$fn > tmpfn; mv tmpfn \$fn; done
  join $a $b -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined1
  join joined1 $c -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined2
  join joined2 $d -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined3
  join joined3 $e -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined4
  join joined4 $f -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined5
  join joined5 $g -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined6
  grep '^Peptide' joined6 > combined
  grep -v '^Peptide' joined6 >> combined
  """
}


process addLociNovelPeptides{
  
  container 'pgpython'
  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  file x from combined_novelpep_output
  
  output:
  val('Finished validating novel peptides') into novelreport
  file 'novel_peptides.txt' into novpeps_finished
  
  """
  python3 /pgpython/group_novpepToLoci.py  --input $x --output novel_peptides.txt --distance 10kb
  """
}


process prepSpectrumAI {

  container 'pgpython'
  
  input:
  set val(peptype), file(x) from variantpsms
  
  output:
  file 'specai_in.txt' into specai_input
  
  """
  python3 /pgpython/label_sub_pos.py --input_psm $x --output specai_in.txt
  """
}


process SpectrumAI {
  container 'spectrumai'

  input:
  file specai_in from specai_input
  file x from specaimzmls

  output: file 'specairesult.txt' into specai

  """
  mkdir mzmls
  for fn in $x; do ln -s `pwd`/\$fn mzmls/; done
  ls mzmls
  Rscript /SpectrumAI/SpectrumAI.R mzmls $specai_in specairesult.txt
  """
}


process SpectrumAIOutParse {

  container 'pgpython'
  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  file x from specai
  file peptides from peptable
  file cosmic
  file dbsnp
  
  output:
  val('Validated variant peptides') into variantreport
  file "variant_peptides.txt" into varpeps_finished
  file "variant_peptides.saav.pep.hg19cor.vcf" into saavvcfs_finished

  """
  python3 /pgpython/parse_spectrumAI_out.py --spectrumAI_out $x --input $peptides --output variant_peptides.txt
  python3 /pgpython/map_cosmic_snp_tohg19.py --input variant_peptides.txt --output variant_peptides.saav.pep.hg19cor.vcf --cosmic_input $cosmic --dbsnp_input $dbsnp
  """
}

variantreport
  .mix(novelreport)
  .subscribe { println(it) }
