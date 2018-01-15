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
repodir = file('.')


/* PIPELINE START */

process prepareContainers {

  """
  workdir=`pwd`
  cd $repodir
  docker build -f $repodir/spectrumAI_Dockerfile -t spectrumai .
  """
}

Channel.fromPath(params.psms).into{novelpsms; variantpsms}


process SpectrumAI {

  container 'spectrumai'
  
  input:
  file x from variantpsms
  
  output:
  file 'parsed_specai.txt'
  
  """
  head -n 1 $x > novelpsms.txt
  egrep '(PGOHUM|lnc)' $x >> novelpsms.txt
  python2.7 label_sub_pos.py --input_psm novelpsms.txt --input_pep example_vardb_6rf_novpep.hg19cor.blastp.annovar.txt --peptide_column "Peptide" --output example_novpep_1mismatch.psm.txt
  RScript 
  python2.7 /Z/jorrit/proteogenomics_python/parse_spectrumAI_out.py --spectrumAI_out specAIout.txt --input example_vardb_6rf_novpep.hg19cor.blastp.annovar.txt --output parsed_specai.txt
  """
}


process novelPsmTable {

  input:
  file x from novelpsms
  
  output:
  file 'novel.fa' into novelfasta
  
  """
  head -n 1 $x > novelpsms.txt
  egrep '(pgohum|lnc)' $x >> novelpsms.txt
  python2.7 /Z/jorrit/proteogenomics_python/to_fasta.py --pep_seq_column 12 --input novelpsms.txt --output novel.fa
  """
}


process BlastPNovel {

  container 'quay.io/biocontainers/blast:2.7.1--boost1.64_1'

  input:
  file novelfasta from novelfasta

  output:
  file 'blastp_out.txt' into novelblast
  
  """
  blastp -db $params.blastdb -query $novelfasta -outfmt '6 qseqid sseqid pident qlen slen qstart qend sstart send mismatch positive gapopen gaps qseq sseq evalue bitscore' -num_threads 8 -max_target_seqs 1 -evalue 1000 -out blastp_out.txt
  """
}

/*
process variantPsmTable {

  input:
  file x from variantpsms
  
  output:
  file 'varpsms.txt' into variantpeptides
  
  """
  echo jellop
  """
}
*/




