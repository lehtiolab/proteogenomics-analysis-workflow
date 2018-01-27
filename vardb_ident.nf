/*
vim: syntax=groovy
-*- mode: groovy;-*-

HiRIEF II varDB pipeline
FIXME:
 - make sure docker files are ok and in right place etc. Try them!
 - nextflow version checker (min 0.26 for transpose to work)
 - merge ident and post cleanly
 - update readme
 - have fasta etc files ready for download
 - default params for everything?
 - check all FIXMEs
 - add authors in NF pipe, all commits count
*/

nf_required_version = '0.26.0'
if( ! nextflow.version.matches(">= ${nf_required_version}") ){
  println("Nextflow version too old, ${nf_required_version} required")
  exit(1)
}


/* SET DEFAULT PARAMS */
mods = file('Mods.txt')
params.ppoolsize = 2
params.gtf = 'vardb.gtf'
params.blastdb = 'Uniprot.Ensembl.RefSeq.GENCODE.proteins.fa'
params.snpfa = 'MSCanProVar_ensemblV79.filtered.fa'
params.genome = 'hg19.fa'

knownproteins = file(params.knownproteins)
blastdb = file(params.blastdb)
gtffile = file(params.gtf)
snpfa = file(params.snpfa)
dbsnp = file(params.dbsnp)
cosmic = file(params.cosmic)
genomefa = file(params.genome)
tdb = file(params.tdb)
ddb = file(params.ddb)
/* PIPELINE START */

Channel
  .from(['target', file(params.tdb)], ['decoy', file(params.ddb)])
  .set { dbs }

Channel
  .fromPath(params.mzmlfile)
  .splitText()
  .map { it -> it.tokenize() }
  .tap { countmzmls; countsets }
  .map { it -> [it[1], it[2], it[0].replaceFirst(/.*fr(\d\d).*/, "\$1").toInteger(), it[0].replaceFirst(/.*\/(\S+)\.mzML/, "\$1"), file(it[0])] }
  .tap { mzmls }
  .combine(dbs)
  .set { dbmzmls }



countmzmls
  .count()
  .set{ amount_mzml }

countsets
  .countBy{ it[1] }
  .set{ amount_sets }

mzmlsets = Channel.create()
mzmlfiles = Channel.create()
mzmls
  .separate(mzmlsets, mzmlfiles) { it -> [it[0], it[4]] }

mzmlfiles
  .collect()
  .into { mzmlfiles_all; specaimzmls; singlemismatch_nov_mzmls }
mzmlsets
  .collect()
  .set { mzmlsets_all }

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
process createSpectraLookup {

  container 'quay.io/biocontainers/msstitch:2.5--py36_0'

  input:
  file mzmlfiles_all
  val mzmlsets_all
  
  output:
  file 'mslookup_db.sqlite' into spec_lookup

  """
  msslookup spectra -i ${mzmlfiles_all.join(' ')} --setnames  ${mzmlsets_all.join(' ')}
  """
}


process msgfPlus {

  /* Latest version has problems when converting to TSV, possible too long identifiers 
     So we use an older version
     LATEST TESTED: container 'quay.io/biocontainers/msgf_plus:2017.07.21--py27_0'
  */
  container 'quay.io/biocontainers/msgf_plus:2016.10.26--py27_1'

  input:
  set val(setname), val(order), val(fraction), val(sample), file(x), val(td), file(db) from dbmzmls

  output:
  set val(setname), val(order), val(fraction), val(sample), file("${sample}.mzid"), val(td) into mzids
  set val(setname), val(td), val(sample), file('out.mzid.tsv') into mzidtsvs
  
  """
  msgf_plus -Xmx16G -d $db -s $x -o "${sample}.mzid" -thread 12 -mod $mods -tda 0 -t 10.0ppm -ti -1,2 -m 0 -inst 3 -e 1 -protocol 4 -ntt 2 -minLength 7 -maxLength 50 -minCharge 2 -maxCharge 6 -n 1 -addFeatures 1
  msgf_plus -Xmx3500M edu.ucsd.msjava.ui.MzIDToTsv -i "${sample}.mzid" -o out.mzid.tsv
  """
}


dmzids = Channel.create()
tmzids = Channel.create()
mzids
  .map { it -> [set:it[0], order:it[1], fr:it[2], sample:it[3], fn:it[4], td:it[5]] }
  .tap { mzids_perco }
  .choice(tmzids, dmzids) { it -> it['td'] == 'target' ? 0 : 1}

tmzids
  .buffer(size: amount_mzml.value)
  .flatMap { it.sort( {a, b -> a['set'] <=> b['set'] ?: a['fr'] <=> b['fr'] }) }
  .buffer { it['order'] == 'last' || it['order'].toInteger() % params.ppoolsize == 0 }
  .map { it -> ["${it[0]['set']}", it.collect() { it['fn'] }, it.collect() { it['sample'] } ] }
  .set { buffer_mzid_target }
dmzids
  .buffer(size: amount_mzml.value)
  .flatMap { it.sort( {a, b -> a['set'] <=> b['set'] ?: a['fr'] <=> b['fr'] }) }
  .buffer { it['order'] == 'last' || it['order'].toInteger() % params.ppoolsize == 0 }
  .map { it -> ["${it[0]['set']}", it.collect() { it['fn'] }, it.collect() { it['sample'] } ] }
  .set { buffer_mzid_decoy }


process percolator {

  container 'quay.io/biocontainers/percolator:3.1--boost_1.623'

  input:
  set val(tset), file('target?'), val(samples) from buffer_mzid_target
  set val(dset), file('decoy?'), val(samples) from buffer_mzid_decoy

  output:
  set val(tset), file('perco.xml') into percolated

  """
  mkdir targets
  mkdir decoys
  count=1;for sam in ${samples.join(' ')}; do ln -s `pwd`/target\$count targets/\${sam}.mzid; echo targets/\${sam}.mzid >> targetmeta; ((count++));done
  count=1;for sam in ${samples.join(' ')}; do ln -s `pwd`/decoy\$count decoys/\${sam}.mzid; echo decoys/\${sam}.mzid >> decoymeta; ((count++));done
  msgf2pin -o percoin.xml -e trypsin targetmeta decoymeta
  percolator -j percoin.xml -X perco.xml --decoy-xml-output -y
  """
}

process filterPercolator {

  container 'quay.io/biocontainers/msstitch:2.5--py36_0'

  input:
  set val(setname), file(x) from percolated
  file 'trypseqdb' from trypseqdb
  file 'protseqdb' from protseqdb
  file knownproteins

  output:
  set val(setname), file('fp_th0.xml') into t_var_filtered_perco
  set val(setname), file('fp_th1.xml') into t_nov_filtered_perco
  set val(setname), file('perco.xml_decoy.xml_h0.xml') into d_nov_filtered_perco
  set val(setname), file('perco.xml_decoy.xml_h1.xml') into d_var_filtered_perco
  """
  msspercolator splittd -i perco.xml 
  msspercolator splitprotein -i perco.xml_target.xml --protheaders '^PGOHUM;^lnc' '^COSMIC;^CanProVar'
  msspercolator splitprotein -i perco.xml_decoy.xml --protheaders '^decoy_PGOHUM;^decoy_lnc' '^decoy_COSMIC;^decoy_CanProVar'
  msspercolator filterseq -i perco.xml_target.xml_h0.xml -o fs_th0.xml --dbfile trypseqdb --insourcefrag 2 --deamidate 
  msspercolator filterseq -i perco.xml_target.xml_h1.xml -o fs_th1.xml --dbfile trypseqdb --insourcefrag 2 --deamidate 
  msspercolator filterprot -i fs_th0.xml -o fp_th0.xml --fasta $knownproteins --dbfile protseqdb --minlen 8 --deamidate --enforce-tryptic
  msspercolator filterprot -i fs_th1.xml -o fp_th1.xml --fasta $knownproteins --dbfile protseqdb --minlen 8 --deamidate --enforce-tryptic
  """
}

/* Group batches per set */ 
t_nov_filtered_perco
  .groupTuple()
  .map { it -> [it[0], 'target', 'novel', it[1]] }
  .set { t_novgrouped_perco }

t_var_filtered_perco
  .groupTuple()
  .map { it -> [it[0], 'target', 'variant', it[1]] }
  .set { t_vargrouped_perco }

d_nov_filtered_perco
  .groupTuple()
  .map { it -> [it[0], 'decoy', 'novel', it[1]] }
  .set { d_novgrouped_perco }

d_var_filtered_perco
  .groupTuple()
  .map { it -> [it[0], 'decoy', 'variant', it[1]] }
  .set { d_vargrouped_perco }

t_novgrouped_perco
  .mix(t_vargrouped_perco, d_novgrouped_perco, d_vargrouped_perco)
  .set { perco_pre_merge }
  

process percolatorMergeBatches {

  container 'quay.io/biocontainers/msstitch:2.5--py36_0'

  input:
  set val(setname), val(td), val(peptype), file('group?') from perco_pre_merge

  output:
  set val(setname), val(td), val(peptype), file('filtered.xml') into perco_merged
  
  """
  msspercolator merge -i group* -o merged.xml
  msspercolator filteruni -i merged.xml -o filtered.xml --score svm
  """
}

perco_t_merged = Channel.create()
perco_d_merged = Channel.create()

perco_merged
  .choice(perco_t_merged, perco_d_merged) { it -> it[1] == 'target' ? 0 : 1}

perco_t_merged
  .groupTuple(by: [0,1,2])
  .buffer(size: amount_sets.value.size() * 2) /* buffer novel and variant of all sets */
  .flatMap { it.sort( {a, b -> a[0] <=> b[0] ?: a[2] <=> b[2] }) }
  .map{ it -> [it[0], it[2], it[3][0]] }
  .set { perco_t_merged_sorted }
perco_d_merged
  .groupTuple(by: [0,1,2])
  .buffer(size: amount_sets.value.size() * 2) /* buffer novel and variant of all sets */
  .flatMap { it.sort( {a, b -> a[0] <=> b[0] ?: a[2] <=> b[2] }) }
  .map{ it -> [it[0], it[2], it[3][0]] }
  .set { perco_d_merged_sorted }

process getQvalityInput {

  container 'quay.io/biocontainers/msstitch:2.5--py36_0'

  input:
  set val(tset), val(peptype), file('target') from perco_t_merged_sorted
  set val(dset), val(peptype), file('decoy') from perco_d_merged_sorted

  output:
  set val(tset), val(peptype), file('tqpsm.txt'), file('dqpsm.txt'), file('tqpep.txt'), file('dqpep.txt'), file('target'), file('decoy') into qvality_input

  """
  msspercolator qvality -i target --decoyfn decoy --feattype psm -o psmqvality.txt || true
  mv target_qvality_input.txt tqpsm.txt
  mv decoy_qvality_input.txt dqpsm.txt
  msspercolator qvality -i target --decoyfn decoy --feattype peptide -o pepqvality.txt || true
  mv target_qvality_input.txt tqpep.txt
  mv decoy_qvality_input.txt dqpep.txt
  """
}

process qvalityMergedBatches {

  container 'quay.io/biocontainers/percolator:3.1--boost_1.623'

  input:
  set val(set), val(peptype), file('tqpsm'), file('dqpsm'), file('tqpep'), file('dqpep'), file('targetperco'), file('decoyperco') from qvality_input
 
  output:
  set val(set), val(peptype), file('qpsm.out'), file('qpep.out'), file('targetperco'), file('decoyperco') into qvality_output
  """
  qvality tqpsm dqpsm -o qpsm.out
  qvality tqpep dqpep -o qpep.out
  """ 
}
  
process recalculatePercolator {

  container 'quay.io/biocontainers/msstitch:2.5--py36_0'

  input:
  set val(set), val(peptype), file('qpsm'), file('qpep'), file('targetperco'), file('decoyperco') from qvality_output

  output:
  set val(set), val(peptype), file('trecalperco.xml'), file('drecalperco.xml') into recal_perco 

  """
  msspercolator reassign -i targetperco --qvality qpsm --feattype psm -o rec_tpsm
  msspercolator reassign -i rec_tpsm --qvality qpep --feattype peptide -o trecalperco.xml
  msspercolator reassign -i decoyperco --qvality qpsm --feattype psm -o rec_dpsm
  msspercolator reassign -i rec_dpsm --qvality qpep --feattype peptide -o drecalperco.xml
  """
}


recal_perco
  .buffer(size: amount_sets.value.size() * 2)
  .flatMap { it.sort( {a, b -> a[1] <=> b[1] ?: a[0] <=> b[0] }) }
  .into { trecalperco; drecalperco }
trecalperco
  .map { it -> [it[0], it[1], it[2]] }
  .set { tpout_perco }
drecalperco
  .map { it -> [it[0], it[1], it[3]] }
  .set { dpout_perco }

mzids_perco
  .map { it -> [it.set, it.td, it.fn] }
  .groupTuple(by: [0,1])
  .buffer(size: amount_sets.value.size() * 2)
  .flatMap { it.sort( {a, b -> a[0] <=> b[0]}) } 
  .into { novmzids; varmzids }
tpercomzids = Channel.create()
dpercomzids = Channel.create()
novmzids
  .concat(varmzids)
  .choice(tpercomzids, dpercomzids) { it -> it[1] == 'target' ? 0 : 1}


process poutToMzidTarget {

  container 'quay.io/biocontainers/pout2mzid:0.3.03--boost1.62_2'

  input:
  set val(pset), val(peptype), file('perco') from tpout_perco
  set val(mset), val(td), file(mzids) from tpercomzids
 
  output:
  set val(pset), val(peptype), file('p2mzid/*.mzid') into tpmzid
  
  """
  ls *.mzid > infiles.txt
  pout2mzid -p perco -i . -f infiles.txt -o p2mzid -c _perco -v
  """
}


process poutToMzidDecoy {

  container 'quay.io/biocontainers/pout2mzid:0.3.03--boost1.62_2'

  input:
  set val(pset), val(peptype), file('perco') from dpout_perco
  set val(mset), val(td), file(mzids) from dpercomzids
 
  output:
  set val(pset), val(peptype), file('p2mzid/*.mzid') into dpmzid
  
  """
  ls *.mzid > infiles.txt
  pout2mzid -p perco -i . -f infiles.txt -o p2mzid -c _perco -v -d
  """
}


varmzidp = Channel.create()
novmzidp = Channel.create()
mzidtsvs
  .buffer(size: amount_mzml.value * 2)
  .flatMap { it.sort( {a, b -> a[1] <=> b[1] ?: a[0] <=> b[0] ?: a[2] <=> b[2]}) }
  .set { sortedtsvs }
tpmzid
  .map { it -> it[2] instanceof List ? it : [it[0], it[1], [it[2]]] }
  .buffer(size: amount_sets.value.size() * 2)
  .flatMap { it.sort( {a, b -> a[0] <=> b[0] ?: a[2] <=> a[2]}) }
  .transpose()
  .set { sorted_tpmzid }
dpmzid
  .map { it -> it[2] instanceof List ? it : [it[0], it[1], [it[2]]] }
  .buffer(size: amount_sets.value.size() * 2)
  .flatMap { it.sort( {a, b -> a[0] <=> b[0] ?: a[2] <=> a[2]}) }
  .transpose()
  .concat(sorted_tpmzid)
  .choice(varmzidp, novmzidp) { it -> it[1] == 'variant' ? 0 : 1}



process annotateMzidTSVPercolator {
  
  container 'quay.io/biocontainers/msstitch:2.5--py36_0'
  
  input:
  set val(vset), val(vartype), file('varmzid') from varmzidp
  set val(nset), val(novtype), file('novmzid') from novmzidp
  set val(pset), val(td), val(sample), file(psms) from sortedtsvs
  
  output:
  set val(td), file("${sample}.txt") into psmsperco

  """
  cp $psms varpsms
  msspsmtable percolator -i $psms -o novperco --mzid novmzid
  msspsmtable percolator -i varpsms -o varperco --mzid varmzid
  cat novperco <( tail -n+2 varperco) > ${sample}.txt 
  """
}

psmsperco
  .filter { it[0] == 'target' }
  .map { it -> it[1] }
  .collect()
  .set { prepsmtable }

process createPSMTable {

  container 'quay.io/biocontainers/msstitch:2.5--py36_0'

  input:
  file 'psms' from prepsmtable
  file 'lookup' from spec_lookup

  output:
  file 'psmtable.txt' into psmtable
  file 'tmp_*' into prepeptables

  """
  msspsmtable merge -o psms.txt -i psms*
  msspsmtable conffilt -i psms.txt -o filtpsm --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'PSM q-value'
  msspsmtable conffilt -i psms.txt -o filtpep --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'peptide q-value'
  cp lookup psmlookup
  msslookup psms -i filtpep --dbfile psmlookup
  msspsmtable specdata -i filtpep --dbfile psmlookup -o psmtable.txt
  sed 's/\\#SpecFile/SpectraFile/' -i psmtable.txt
  mkdir splitpsms
  msspsmtable split -i psmtable.txt --splitcol 2 -d splitpsms
  for fn in `ls splitpsms`;do msspeptable psm2pep -i splitpsms/\$fn -o tmp_\$fn --scorecolpattern svm --spectracol 1; done
  """
}

prepeptables
  .flatten()
  .set { prepeptable }

process peptable {

  container 'ubuntu:latest'

  input:
  file x from prepeptable

  output:
  file "${x.baseName.replaceFirst(/tsv/, '')}.txt" into peptable

  """
  paste <( cut -f 12 $x) <( cut -f 1-11,13-22 $x) > ${x.baseName.replaceFirst(/tsv/, '')}.txt
  """
}

peptable
  .collect()
  .set { peptables }


process SplitPSMTableNovelVariant {
  input:
  file x from psmtable
  
  output:
  file 'variantpsms' into variantpsms
  file 'novelpsms' into novelpsms

  """
  head -n 1 $x > variantpsms
  head -n 1 $x > novelpsms
  egrep '(PGOHUM|lnc)' $x >> novelpsms
  egrep '(COSMIC|CanProVar)' $x >> variantpsms
  """
}

novelpsms
  .into{novelpsmsFastaBedGFF; novelpsms_specai}


process createFastaBedGFF {
 container 'pgpython'

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


/* FIXME this needs to be made conditional */

bamFiles = Channel
  .fromPath(params.bamfiles)
  .map { fn -> [ fn, fn + '.bai' ] }
  .collect()

process scanBams {
  container 'pgpython'

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

/*
process makePlots {

}

*/

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

process combineResults{
  
  container 'pgpython'

  input:
  file a from ns_snp_out
  file b from novpep_singlemisspecai
  file c from peptable_blat
  file d from annovar_parsed
  file e from phastcons_out
  file f from phylocsf_out
  file g from scannedbams
  
  output:
  file 'combined' into combined_novelpep_output
  
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
  file x from variantpsms
  
  output:
  file 'specai_in.txt' into specai_input
  
  """
  head -n 1 $x > variantpsms.txt
  egrep '(COSMIC|CanProVar)' $x >> variantpsms.txt
  python3 /pgpython/label_sub_pos.py --input_psm variantpsms.txt --output specai_in.txt
  """
}


process SpectrumAI {
  container 'spectrumai'

  input:
  file specai_in from specai_input
  file x from specaimzmls

  output: file 'specairesult.txt' into specai

  """
  #mkdir mzmls
  #for fn in $x; do ln -s `pwd`/\$fn mzmls/; done
  mkdir mzmls
  cd mzmls
  for fn in $x; do ln -s ../\$fn .; done
  cd ..
  ls mzmls
  Rscript /SpectrumAI/SpectrumAI.R mzmls $specai_in specairesult.txt
  """
}


process SpectrumAIOutParse {

  container 'pgpython'
  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  file x from specai
  each peptides from peptables
  file cosmic
  file dbsnp
  
  output:
  val('Validated variant peptides of a sample set') into variantreport
  file "variant_peptides_${peptides.baseName}" into varpeps_finished
  file "saav.pep.hg19cor.${peptides.baseName.replaceFirst(/\.txt$/, '')}.vcf" into saavvcfs_finished

  """
  python3 /pgpython/parse_spectrumAI_out.py --spectrumAI_out $x --input $peptides --output variant_peptides_${peptides.baseName}
  python3 /pgpython/map_cosmic_snp_tohg19.py --input variant_peptides_${peptides.baseName} --output saav.pep.hg19cor.${peptides.baseName.replaceFirst(/\.txt$/, '')}.vcf --cosmic_input $cosmic --dbsnp_input $dbsnp
  """
}

variantreport
  .mix(novelreport)
  .subscribe { println(it) }
