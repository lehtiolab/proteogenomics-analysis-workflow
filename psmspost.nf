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
params.genomeFasta = 'hg19.fa'
gtffile = file(params.gtf)
fafile = file(params.fasta)

repodir = file('.')


/* PIPELINE START */

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

Channel.fromPath(params.psms).into{novelpsms; variantpsms}
containers_done = Channel.from(1)

process createFastaBedGFF {
 container 'pgpython'

 input:
 file x from novelpsms
 val tf from containers_done

 output:
 file 'novel_peptides.fa' into novelfasta
 file 'novel_peptides.bed' into novelbed
 file 'novel_peptides.gff3' into novelGFF3
 file 'novel_peptides.tab.txt' into novelpep

 """
 head -n 1 $x > novelpsms.txt
 egrep '(PGOHUM|lnc)' $x >> novelpsms.txt
 python3 /pgpython/map_novelpeptide2genome.py --input novelpsms.txt --gtf $gtffile --fastadb $fafile --tab_out novel_peptides.tab.txt --fasta_out novel_peptides.fa --gff3_out novel_peptides.gff3 --bed_out novel_peptides.bed

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

process ParseBlastpOut {
 container 'pypython'
 
 input:
 file novelpep from novelpep
 file novelblast from novelblast

 output:
 file 'peptable_blastp' into peptableBlastp

 """
 python parse_blastp_output.py --input $novelpep --blastp_result novelblast --fasta $params.blastdb --output peptable_blastp.txt

 """

}

/*
process BLATNovel {
  container ''

  input:
  file novelfasta from novelfasta

  output:
  file 'blat_out.pslx' into novelblat

  """
  blat $params.genomeFasta $novelfasta -t=dnax -q=prot -tileSize=5 -minIdentity=99 -out=pslx blat_out.pslx 
  """
}

process parseBLATout {
 container 'pypython'

 input:
 file novelblat from novelblat
 file novelpep from novelpep

 output:
 file 'peptable_blat.txt' into peptableBlat

 """
 python parse_BLAT_out.py $novelblat $novelpep peptable_blat.txt

 """
}

process labelnsSNP {

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

process annovar {
  
  container 'annovar'
  
  input:
  file novelbed
  output:
  file 'novpep_annovar.variant_function' into annovar

  """
  /annovar/annotate_variation.pl -out novpep_annovar -build hg19 $novelbed /annovar/humandb/
  """

}

