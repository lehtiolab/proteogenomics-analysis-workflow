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

Channel.fromPath(params.mzmls).into{tmzmls; dmzmls}

process msgfPlusTarget {

  container 'quay.io/biocontainers/msgf_plus:2017.07.21--py27_0'

  input:
  file x from tmzmls

  output:
  file 'toutmzid.mzid' into tmzids
  file 'touttsv.tsv' into tmzidtsvs
  
  """
  msgf_plus -Xmx16G -d $tdb -s $x -o toutmzid.mzid -thread 12 -mod $mods -tda 0 -t 10.0ppm -ti -1,2 -m 0 -inst 3 -e 1 -protocol 4 -ntt 2 -minLength 7 -maxLength 50 -minCharge 2 -maxCharge 6 -n 1 -addFeatures 1
  msgf_plus -Xmx3500M edu.ucsd.msjava.ui.MzIDToTsv -i toutmzid.mzid -o touttsv.tsv
  """
}


process msgfPlusDecoy {

  container 'quay.io/biocontainers/msgf_plus:2017.07.21--py27_0'

  input:
  file x from dmzmls

  output:
  file 'doutmzid.mzid' into dmzids
  file 'douttsv.tsv' into dmzidtsvs
  
  """
  msgf_plus -Xmx16G -d $ddb -s $x -o doutmzid.mzid -thread 12 -mod $mods -tda 0 -t 10.0ppm -ti -1,2 -m 0 -inst 3 -e 1 -protocol 4 -ntt 2 -minLength 7 -maxLength 50 -minCharge 2 -maxCharge 6 -n 1 -addFeatures 1
  msgf_plus -Xmx3500M edu.ucsd.msjava.ui.MzIDToTsv -i doutmzid.mzid -o douttsv.tsv
  """
}


bufferedtmzids = tmzids.tap{ singletmzids }.buffer(size: params.ppoolsize, remainder: true)
buffereddmzids = dmzids.tap{ singledmzids }.buffer(size: params.ppoolsize, remainder: true)


process percolator {
  
  /* container 'quay.io/biocontainers/percolator:3.1--boost_1.622' */
  container 'quay.io/biocontainers/percolator:3.1--boost_1.622'
  

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

