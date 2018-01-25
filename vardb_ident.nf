/*
vim: syntax=groovy
-*- mode: groovy;-*-

HiRIEF II varDB pipeline
post identification and FDR

incorporate to main pipe later

*/


/* SET DEFAULT PARAMS */
mods = file('Mods.txt')
params.ppoolsize = 2
knownproteins = file(params.knownproteins)

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
  .set{ dbmzmls }

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
  .set { mzmlfiles_all }
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
  ls
  count=1;for sam in ${samples.join(' ')}; do ln -s `pwd`/target\$count targets/\${sam}.mzid; echo targets/\${sam}.mzid >> targetmeta; ((count++));done
  count=1;for sam in ${samples.join(' ')}; do ln -s `pwd`/target\$count decoys/\${sam}.mzid; echo decoys/\${sam}.mzid >> decoymeta; ((count++));done
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
  file "${sample}.txt" into psmsperco
  """
  cp $psms varpsms
  msspsmtable percolator -i $psms -o novperco --mzid varmzid
  msspsmtable percolator -i varpsms -o varperco --mzid novmzid
  cat novperco <( tail -n+2 varperco) > ${sample}.txt 
  """
}

psmsperco
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
  file 'peptable.txt' into peptable

  """
  paste <( cut -f 12 $x) <( cut -f 1-11,13-22 $x) > peptable.txt
  """
}
