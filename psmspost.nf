/*
vim: syntax=groovy
-*- mode: groovy;-*-

HiRIEF II varDB pipeline
post identification and FDR

incorporate to main pipe later
dependencies in special dir

*/


/* SET DEFAULT PARAMS */
params.blastdb = 'UniProteome+Ensembl87+refseq+GENCODE24.proteins.fasta'
params.variantout = 'variant_peptides.txt'
gtffile = file(params.gtf)
fafile = file(params.fasta)

repodir = file('.')


/* PIPELINE START */
specaimzmls = Channel.fromPath(params.mzmls).collect()
peptidetable = Channel.fromPath(params.peptable)

/*
process prepareContainers {

  output: val 1 into containers_done
  
  """
  workdir=`pwd`
  cd $repodir
  docker build -f $repodir/spectrumAI_Dockerfile -t spectrumai .
  docker build -f $repodir/pgpython_Dockerfile -t pgpython .
  """
}
*/
containers_done = Channel.from(1)
Channel.fromPath(params.psms).into{novelpsms; variantpsms}


process prepSpectrumAI {

  container 'pgpython'
  
  input:
  val tf from containers_done
  file x from variantpsms
  
  output:
  file 'specai_in.txt' into specai_input
  
  """
  head -n 1 $x > variantpsms.txt
  egrep '(COSMIC|CanProVar)' $x >> variantpsms.txt
  python2.7 /pgpython/label_sub_pos.py --input_psm variantpsms.txt --output specai_in.txt
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
  cd mzmls
  for fn in $x; do ln -s ../\$fn .; done
  cd ..
  ls mzmls
  Rscript /SpectrumAI/SpectrumAI.R mzmls $specai_in specairesult.txt
  """
}


process SpectrumAIOutParse {

  container 'pgpython'

  input:
  file x from specai
  file 'peptide_table.txt' from peptidetable
  
  output:
  file 'parsed_specai.txt' into variantpep_output

  """
  python2.7 /pgpython/parse_spectrumAI_out.py --spectrumAI_out $x --input peptide_table.txt --output parsed_specai.txt
  """
}


variantpep_output
  .collectFile(name: file(params.variantout))
  .println {"Variant peptides saved to file: $it" }
 
