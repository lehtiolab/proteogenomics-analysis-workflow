/*
vim: syntax=groovy
-*- mode: groovy;-*-

HiRIEF II varDB pipeline
post identification and FDR

incorporate to main pipe later
dependencies in special dir

*/


/* SET DEFAULT PARAMS */
params.genomeFasta = 'hg19.fa'

blastdb = file(params.blastdb)
gtffile = file(params.gtf)
snpdb = file(params.snpdb)
genomefa = file(params.genomeFasta)
fafile = file(params.fasta)

repodir = file('.')


/* PIPELINE START */
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

Channel.fromPath(params.psms).into{novelpsms; variantpsms}
containers_done = Channel.from(1)

process createFastaBedGFF {
 container 'pgpython'

 input:
 file x from novelpsms
 file gtffile
 file fafile
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
 file novelpep from blastnovelpep
 file novelblast from novelblast
 file blastdb

 output:
 file 'peptable_blastp.txt' into peptable_blastp

 """
 python3 /pgpython/parse_BLASTP_out.py --input $novelpep --blastp_result $novelblast --fasta $blastdb --output peptable_blastp.txt
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
  file snpdb

  output:
  file 'nssnp.txt' into ns_snp_out

  """
  echo hello
  python3 /pgpython/label_nsSNP_pep.py --input $peptable --nsSNPdb $snpdb --output nssnp.txt
  """
}

novelGFF3
  .into { novelGFF3_phast; novelGFF3_phylo; novelGFF3_bams }

process phastcons {
  container 'phastcons'
  
  input:
  file novelgff from novelGFF3_phast
  output:
  file 'phastcons.txt' into phastcons_out

  """
  python3 /pgpython/calculate_phastcons.py $novelgff /hg19.100way.phastCons.bw phastcons.txt
  """
}

process phyloCSF {
  
  container 'phylocsf'

  input:
  file novelgff from novelGFF3_phylo

  """
  python3 /pgpython/calculate_phylocsf.py $novelgff > phylocsf.txt
  """

}
/*
process scanBams {
  container 'pgpython'

  input:
  file gff from novelGFF3_bams
  file bams from bamFiles
  
  output:
  file 'scannedbams.txt' into scannedbams

  """
  python3 /pgpython/scan_bams.py  --gff_input $gff --bam_files $bams --output scannedbams.txt
  """
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
