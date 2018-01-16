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

  """
  echo ha
  python2.7 /pgpython/parse_spectrumAI_out.py --spectrumAI_out $x --input peptide_table.txt --output parsed_specai.txt
  echo he
  """
}


/*

process createFasta.Bed.GFF.txt {
 container 'pgpython'
 input:
 file x from novelpsms
 val tf from containers_done

 output:
 file 'novel_peptides.fa' into novelfasta
 file 'novel_peptides.bed' into novelbed
 file 'novel_peptides.gff3' into novelGFF3
 file 'novep_peptides.tab.txt' into novelpep

 """
 head -n 1 $x > novelpsms.txt
 egrep '(PGOHUM|lnc)' $x >> novelpsms.txt

 python2.7 /pgpython/map_novelpeptide2genome.py --input $x --gtf $gtffile --fastadb $fafile --tab_out novel_peptides.tab.txt --fasta_out novel_peptides_fa --gff3_out novel_peptides.gff3 --bed_out novel_peptides.bed

 """
}
*/

/*
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

process BLATNovel {
  container ''

  input:
  file novelfasta from novelfasta

  output:
  file 'blat_out.txt' into novelblat

  """
  blat -db $params.genomedb -query $novelfasta 

}

process labelnsSNP {

}

process anovar {

}

process phastcon {

}

process phyloCSF {

}

process scanBam {

}

process parseBlastpOutput{

}

process parseBlatOutput{

}

process combineResults{
}

process makePlots {

}

*/
