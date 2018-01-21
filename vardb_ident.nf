/*
vim: syntax=groovy
-*- mode: groovy;-*-

HiRIEF II varDB pipeline
post identification and FDR

incorporate to main pipe later

*/


/* SET DEFAULT PARAMS */
mods = file('Mods.txt')
params.ppoolsize = 8
tdb = file(params.tdb)
ddb = file(params.ddb)


/* PIPELINE START */

Channel
  .fromPath(params.mzmls)
  .map { it -> ['setA', it.name.replace('.mzML', ''), it] }
  .into{ dmzmls; tmzmls }

process Test {
  input:
  set val(setname), val(filename), file(x) from tmzmls
  output:
  set val(setname), val(filename), file(x) into testfns
  """
  echo hei
  """
}

testfns
  .map{ it -> [set: it[0], sample: it[1], file: it[2]] }
  .buffer(size: 6)
  .flatMap{ it.sort( {a, b -> a['set'] <=> b['set'] ?: a['sample'] <=> b['sample'] }) }
  .buffer(size: params.ppoolsize, remainder: true)
  .subscribe{ println(it) }

/*
process msgfPlusTarget {

  container 'quay.io/biocontainers/msgf_plus:2017.07.21--py27_0'

  input:
  set val(setname), file(x) from tmzmls
  file tdb

  output:
  file 'toutmzid.mzid' into tmzids
  file 'touttsv.tsv' into tmzidtsvs
  
  """
  msgf_plus -Xmx16G -d $tdb -s $x -o toutmzid.mzid -thread 12 -mod $mods -tda 0 -t 10.0ppm -ti -1,2 -m 0 -inst 3 -e 1 -protocol 4 -ntt 2 -minLength 7 -maxLength 50 -minCharge 2 -maxCharge 6 -n 1 -addFeatures 1
  msgf_plus -Xmx3500M edu.ucsd.msjava.ui.MzIDToTsv -i toutmzid.mzid -o touttsv.tsv
  """
}
*/
process msgfPlusDecoy {

  container 'quay.io/biocontainers/msgf_plus:2017.07.21--py27_0'

  input:
  set val(setname), val(sample), file(x) from dmzmls
  file ddb

  output:
  set val(setname), val(sample), file('doutmzid.mzid') into dmzids
  set val(setname), val(sample), file('douttsv.tsv') into dmzidtsvs
  
  """
  msgf_plus -Xmx16G -d $ddb -s $x -o doutmzid.mzid -thread 12 -mod $mods -tda 0 -t 10.0ppm -ti -1,2 -m 0 -inst 3 -e 1 -protocol 4 -ntt 2 -minLength 7 -maxLength 50 -minCharge 2 -maxCharge 6 -n 1 -addFeatures 1
  msgf_plus -Xmx3500M edu.ucsd.msjava.ui.MzIDToTsv -i doutmzid.mzid -o douttsv.tsv
  """
}

dmzids
  .map{ it -> [set: it[0], sample: it[1], file: it[2]] }
  .buffer(size: 6)
  .flatMap{ it.sort( {a, b -> a['set'] <=> b['set'] ?: a['sample'] <=> b['sample'] }) }
  .tap { singledmzids }
  .buffer(size: params.ppoolsize, remainder: true)
  .subscribe{ println(it) }

/*

process percolator {
  
  container 'quay.io/biocontainers/percolator:3.1--boost_1.623'

  input:
  file 'tmzid' from bufferedtmzids
  file 'dmzid' from buffereddmzids
  
  output:
  file 'perco.xml' into percolateds
  
  """
  echo tmzid* | sed 's/ /\\n/g' > targetmeta
  echo dmzid* | sed 's/ /\\n/g' > decoymeta
  export LC_ALL=en_US.UTF-8 
  msgf2pin -o percoin.xml -e trypsin targetmeta decoymeta
  percolator -j percoin.xml -X perco.xml --decoy-xml-output -y
  """
}
*/


/*
process collectPSMFractions {
  input:
  file 'psms_*' from tmzidtsvs.collect()

  output:
  file 'mergedpsms.txt' into mergedpsms

  """
  echo psms_* && head -n 1 psms_1 > mergedpsms.txt && tail -n+2 psms_* >> mergedpsms.txt
  """
}
process createNovelPeptideTable {
  input:
  file psms from mergedpsms

  output:
  file 'peptides.txt' into novelpeptides
  

  """
  source $params.venv/bin/activate
  msspeptable psm2pep -i $psms -o peptides.txt --spectracol 1 --scorecolpattern 'MSGFScore' 
  """
}

*/

