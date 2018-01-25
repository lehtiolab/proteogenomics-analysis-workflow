/*
vim: syntax=groovy
-*- mode: groovy;-*-

HiRIEF II varDB pipeline
post identification and FDR

incorporate to main pipe later
dependencies in special dir

*/


/* SET DEFAULT PARAMS */
params.gtf = 'vardb.gtf'
params.blastdb = 'Uniprot.Ensembl.RefSeq.GENCODE.proteins.fa'
params.snpfa = 'MSCanProVar_ensemblV79.filtered.fa'
params.genome = 'hg19.fa'


blastdb = file(params.blastdb)
gtffile = file(params.gtf)
snpfa = file(params.snpfa)
dbsnp = file(params.dbsnp)
cosmic = file(params.cosmic)
genomefa = file(params.genome)
fafile = file(params.fasta)

repodir = file('.')


/* PIPELINE START */

psms = Channel.fromPath(params.psms)
peptidetable = Channel.fromPath(params.peptable)
Channel
  .fromPath(params.mzmls)
  .collect()
  .into { specaimzmls; singlemismatch_nov_mzmls }


process SplitPSMTableNovelVariant {
  input:
  file x from psms
  
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
 file fafile

 output:
 file 'novel_peptides.fa' into novelfasta
 file 'novel_peptides.bed' into novelbed
 file 'novel_peptides.gff3' into novelGFF3
 file 'novel_peptides.tab.txt' into novelpep

 """
 python3 /pgpython/map_novelpeptide2genome.py --input $novelpsmsFastaBedGFF --gtf $gtffile --fastadb $fafile --tab_out novel_peptides.tab.txt --fasta_out novel_peptides.fa --gff3_out novel_peptides.gff3 --bed_out novel_peptides.bed
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
  cd mzmls
  for fn in $mzml; do ln -s ../\$fn .; done
  cd ..
  Rscript /SpectrumAI/SpectrumAI.R mzmls $x singlemis_specai.txt
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
  file 'peptide_table.txt' from peptidetable
  file cosmic
  file dbsnp
  
  output:
  val('Finished validating variant peptides') into variantreport
  file 'variant_peptides.txt' into varpeps_finished
  file 'saav.pep.hg19cor.txt' into saavvcfs_finished

  """
  python3 /pgpython/parse_spectrumAI_out.py --spectrumAI_out $x --input peptide_table.txt --output variant_peptides.txt
  python3 /pgpython/map_cosmic_snp_tohg19.py --input variant_peptides.txt --output saav.pep.hg19cor.txt --cosmic_input $cosmic --dbsnp_input $dbsnp
  """
}

variantreport
  .mix(novelreport)
  .subscribe { println(it) }
