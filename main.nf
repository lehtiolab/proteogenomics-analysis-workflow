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



/* PIPELINE START */

// Either feed an mzmldef file (tab separated lines with filepath\tsetname), or /path/to/\*.mzML
if (!params.mzmldef && !params.input) {
  Channel
    .fromPath(params.mzmls)
    .map { it -> [it, 'NA'] }
    .set { mzml_in }
} else {
  header = ['mzmlfile', 'setname', 'plate', 'fraction']
  mzmldef = params.mzmldef ?: params.input
  mzmllines = file(mzmldef).readLines().collect { it.tokenize('\t') }
  if (mzmllines[0] == header) {
    /* As above, future use with pushing files with a header becomes enabled, as long as
    they use this header format. We cannot do module importing etc yet, have to use DSL2
    for that. That is something to strive for in the future.
    */
    mzmllines.remove(0)
  }
  Channel
    .from(mzmllines)
    .set { mzml_in }
}


// Isobaric input parsing to setisobaric and setdenoms maps
// example: --isobaric 'set1:tmt10plex:127N:128N set2:tmtpro:sweep set3:itraq8plex:intensity'
isop = params.isobaric ? params.isobaric.tokenize(' ') : false
setisobaric = isop ? isop.collect() {
  y -> y.tokenize(':')
}.collectEntries() {
  x-> [x[0], x[1]]
} : false
// FIXME add non-isobaric sets here if we have any mixed-in?
setdenoms = isop ? isop.collect() {
  y -> y.tokenize(':')
}.collectEntries() {
  x-> [x[0], x[2..-1]]
} : false


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
  msstitch split -i normalpsms --splitcol bioset
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
  msstitch split -i $normpsm --splitcol `python -c 'with open("$normpsm") as fp: h=next(fp).strip().split("\\t");print(h.index("Strip")+1)'`
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
  msstitch peptides -i $normpsm -o peptides --scorecolpattern area --spectracol 1 
  """
}

pipep = params.pisepdb ? Channel.fromPath(params.pisepdb) : Channel.empty()
varnov_peptides = params.pisepdb ? Channel.fromPath(params.pisepdb) : Channel.from(tdb)
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
  echo \'${groovy.json.JsonOutput.toJson(strip)}\' >> strip.json
  pi_database_splitter.py -i $pipeptides -p $peptides --stripdef strip.json --deltacolpattern Delta --fraccolpattern Fraction --fdrcolpattern '^q-value' --picutoff 0.2 --fdrcutoff 0.0 --maxlen $params.maxlen --minlen $params.minlen
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

process makeTargetSeqLookup {

  input:
  file(tdb) from varnov_peptides
  file(knownproteins)

  output:
  set file('mslookup_db.sqlite'), file('decoy_known.fa') into target_seq_lookup

  script:
  """
  # create text file of pi sep (much faster to import to SQLite than fasta)
  ${params.pisepdb ? "cut -f2 $tdb | sort -u > targetseq.txt" : "grep -v '^>' $tdb > targetseq.txt"}
  # Add trypsinized known proteins to txt file
  msstitch trypsinize -i $knownproteins -o knowntryp
  grep -v '^>' knowntryp >> targetseq.txt

  # TODO parametrize notrypsin?
  msstitch storeseq -i targetseq.txt --minlen $params.minlen ${params.pisepdb ? '--notrypsin': ''}
  msstitch makedecoy -i $knownproteins --dbfile mslookup_db.sqlite -o decoy_known.fa --scramble tryp_rev --minlen $params.minlen
  """
}

// Join DB with mzmls, use join filters out DBs without a match (in case of missing mzML fractions) 
// or duplicates (when having reruns, or when not running pI separared DBs in which case you only need 
// this process once). So we dont generate more decoys than necessary
db_w_id
  .join(mzml_dbfilter)
  .map { it -> [it[0], it[1]] }
  .combine(target_seq_lookup)
  .set { db_filtered }

process concatFasta {
 
  input:
  set val(dbid), file(db), file(targetlookup), file('decoy_known.fa') from db_filtered
  file knownproteins

  output:
  set val(dbid), file('td_concat.fa') into db_concatdecoy

  script:
  """
  # copy DB for faster access on network FS
  cp ${targetlookup} localdb.sql
  cat $db $knownproteins > td_concat.fa
  msstitch makedecoy -i $db --dbfile localdb.sql -o decoy_db.fa --scramble tryp_rev --minlen $params.minlen ${params.pisepdb ? '--notrypsin': ''}
  cat decoy_db.fa decoy_known.fa >> td_concat.fa
  rm decoy_db.fa localdb.sql
  """
}

// Now re-match the DB (now with decoy) with mzML, this is needed to fan out if more mzMLs than
// DBs have been used so we use the cross operator
db_concatdecoy
  .cross(mzml_dbid) // gives two Arrays so unfold them in next map step
  .map { it -> [it[0][0], it[0][1], it[1][1], it[1][2], it[1][3]] } // dbid, db, set, sample, file
  .set { mzml_msgf }

process prepareFilterDB {

  input:
  file(knownproteins)
  file(snpfa)

  output:
  file('knownprot.sqlite') into protseqdb
  file('snpprot.sqlite') into snpseqdb
  file('mslookup_db.sqlite') into trypseqdb

  """
  msstitch storeseq -i $knownproteins --minlen $params.minlen --fullprotein --minlen 7 -o knownprot.sqlite
  msstitch storeseq -i $snpfa --minlen $params.minlen -o snpprot.sqlite --fullprotein --minlen 7
  msstitch storeseq -i $knownproteins --insourcefrag --minlen $params.minlen
  """
}


process IsobaricQuant {

  when: !params.quantlookup && params.isobaric

  input:
  set val(setname), val(sample), file(infile), val(strip), val(fraction) from mzml_isobaric

  output:
  set val(sample), file("${infile}.consensusXML") into isobaricxml

  script:
  activationtype = [hcd:'High-energy collision-induced dissociation', cid:'Collision-induced dissociation', etd:'Electron transfer dissociation'][params.activation]
  isobtype = setisobaric && setisobaric[setname] ? setisobaric[setname] : false
  isobtype = isobtype == 'tmtpro' ? 'tmt16plex' : isobtype
  plextype = isobtype ? isobtype.replaceFirst(/[0-9]+plex/, "") : 'false'
  massshift = [tmt:0.0013, itraq:0.00125, false:0][plextype]
  """
  IsobaricAnalyzer  -type $isobtype -in $infile -out "${infile}.consensusXML" -extraction:select_activation "$activationtype" -extraction:reporter_mass_shift $massshift -extraction:min_precursor_intensity 1.0 -extraction:keep_unannotated_precursor true -quantification:isotope_correction true 
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


process createNewSpectraLookup {

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
  msstitch storespectra --spectra ${mzmlfiles.join(' ')} --setnames ${setnames.join(' ')}
  msstitch storequant --dbfile mslookup_db.sqlite --isobaric ${isobxmls.join(' ')} --spectra ${mzmlfiles.join(' ')}
  """
  else
  """
  msstitch storespectra --spectra ${mzmlfiles.join(' ')} --setnames ${setnames.join(' ')}
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
  set val(setfr_id), file(db), val(setname), val(sample), file(mzml) from mzml_msgf
  file mods

  output:
  set val(setname), val(sample), file("${sample}.mzid") into mzids
  set val(setname), file("${sample}.mzid"), file('out.mzid.tsv') into mzidtsvs
  
  script:
  mem = db.size() * 16 // used in conf profile
  msgfprotocol = 0
  """
  msgf_plus -Xmx${task.memory.toMega()}M -d $db -s $mzml -o "${sample}.mzid" -thread ${task.cpus * params.threadspercore} -mod $mods -maxMissedCleavages ${params.maxmiscleav} -tda 0 -t 10.0ppm -ti -1,2 -m 0 -inst 3 -e 1 -protocol ${msgfprotocol} -ntt 2 -minLength $params.minlen -maxLength $params.maxlen -minCharge 2 -maxCharge 6 -n 1 -addFeatures 1
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

  script:
  mzmlcount = samples.size() 
  """
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

  script:
  """
  msstitch splitperco -i $x --protheaders "known:${params.knownheaders}|novel:${params.varheaders}"
  """
}


process getNovelPercolator {

  when: params.novheaders

  input:
  set val(setname), file(x) from nov_percolated

  output:
  set val(setname), val('novel'), file("${x}_h0.xml") into nov_perco

  script:
  """
  msstitch splitperco -i $x --protheaders "known:${params.knownheaders}|novel:${params.novheaders}"
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
  msstitch filterperco -i $perco -o filtseq --dbfile trypseqdb --insourcefrag 2 --deamidate 
  msstitch filterperco -i filtseq -o filtprot --fullprotein --dbfile protseqdb --minlen $params.minlen --deamidate
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
  set val(setname), file('mzid?'), file('tsv?'), val(peptype), file(perco) from allmzperco 

  output:
  set val(setname), val(peptype), file('target.tsv') into mzidtsv_perco

  script:
  """
  tsvs=""
  mzids=""
  count=1; for tsvfn in \$(ls tsv*)
    do 
    tsvs="\${tsvs} tsv\${count}"
    mzids="\${mzids} mzid\${count}"
    ((count++))
    done
  mkdir outtables
  msstitch perco2psm --perco $perco -d outtables -i \$tsvs --mzids \$mzids
  msstitch concat -i outtables/* -o psms
  msstitch split -i psms --splitcol TD
  """
}

mzidtsv_perco
  .combine(spec_lookup)
  .set { prepsm }

process createPSMTable {

  input:
  set val(setname), val(peptype), file('psms'), file('lookup') from prepsm
  val(mzmlcount) from mzmlcount_psm

  output:
  set val(setname), val(peptype), file("${setname}_${peptype}_psmtable.txt") into psmtable

  script:
  """
  msstitch conffilt -i psms -o filtpsm --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'PSM q-value'
  msstitch conffilt -i filtpsm -o filtpep --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'peptide q-value'
  cp lookup psmlookup
  msstitch psmtable -i filtpep --dbfile psmlookup --addbioset -o ${setname}_${peptype}_psmtable.txt ${params.isobaric ? '--isobaric': ''}
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
  msstitch peptides -i psms.txt -o peptidetable.txt --scorecolpattern svm --spectracol 1 \
    ${setisobaric && setisobaric[setname] ? "--isobquantcolpattern plex --minint 0.1 --logisoquant --denompatterns ${setdenoms[setname].join(' ')}" : ''}
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
  .set { grouped_saavnov_mzml_psms }
  

process ValidateSingleMismatchNovpeps {
  
  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { it == "precursorError.histogram.plot.pdf" ? "${setname}_novel_precursorError_plot.pdf" : it }
  
  input:
  set val(setname), val(samples), file(mzmls), file(psms) from grouped_saavnov_mzml_psms

  output:
  set val(setname), file("${setname}_novel_saav_specai.txt") into singlemis_specai
  file 'precursorError.histogram.plot.pdf' optional true into novel_specai_plot

  """
  mkdir mzmls
  for fn in $mzmls; do ln -s `pwd`/\$fn mzmls/; done
  Rscript /SpectrumAI/SpectrumAI.R mzmls $psms ${setname}_novel_saav_specai.txt || cp $psms singlemis_specai.txt
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
  file(snpdb) from snpseqdb

  output:
  set val(setname), file('nssnp.txt') into ns_snp_out

  """
  label_nsSNP_pep.py --input $peptable --nsSNPdb $snpfa --dbfile "$snpdb" --output nssnp.txt --minlen $params.minlen
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
  fixfields=`head -n1 peps1 |tr -s '\\t' '\\n' | egrep -vn '(Setname|Spectrum|Files|Charge|q-val|plex|${acc_removemap[peptype]})' | cut -f 1 -d ':'`
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
