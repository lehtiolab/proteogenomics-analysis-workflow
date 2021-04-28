Integrated proteogenomics analysis workflow
==============

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/glormph/ipaw.svg)](https://hub.docker.com/r/glormph/ipaw)

This is a workflow to identify, curate, and validate variant and novel peptides from MS proteomics spectra data, using databases containing novel and variant peptides, such as the VarDB database. VarDB combines entries from COSMIC, PGOHUM, CanProVar and lncipedia. The workflow takes mzML spectra files as input, is powered by [Nextflow](https://nextflow.io) and runs in [Docker](https://docker.com) or [Singularity](https://sylabs.io/singularity) containers.

Searches are run using [MSGF+](https://omics.pnl.gov/software/ms-gf) on a concatenated target and decoy databases which are then passed to [Percolator](http://percolator.ms) for statistical evaluation, in which FDR is determined in a class specific manner, filtering out known peptides and dividing novel/variant in different FDR arms. Thereafter a curation procedure is performed in which resulting peptides are evaluated on several different criteria, dependent on the peptide.

Please cite the following paper when you have used the workflow for publications :)

Zhu Y, Orre LM, Johansson HJ, Huss M, Boekel J, Vesterlund M, Fernandez-Woodbridge A, Branca RMM, Lehtio J: Discovery of coding regions in the human genome by integrated proteogenomics analysis workflow. Nat Commun 2018, 9(1):903.  [PMID: 29500430](https://www.ncbi.nlm.nih.gov/pubmed/29500430)

![workflow image](https://github.com/lehtiolab/proteogenomics-analysis-workflow/blob/master/assets/workflow.png)

### Before running 

  + Install [Docker](https://docker.io) or [Singularity](https://sylabs.io/singularity)
  + Install [Nextflow](https://nextflow.io)


### Detailed pipeline inputs

  + Database search related inputs for MSGFplus
    + Spectra files input
    
    `--mzmldef  # a tab deliminated text file with mzmlfilepath(absolute path) and setname`
 
    + Modification file for MSGF+. Default file is for TMT labelled samples. [Here is an example.](https://bix-lab.ucsd.edu/download/attachments/13533355/Mods.txt?version=2&modificationDate=1358975546000)
    
    `--mods Mods.txt # use standard Unimod name for modification`
    
    + Fragment method
    
    `--activation hcd  # default, else use cid, etd`
    
    + Specify search DB
    
    `--tdb /path/to/vardb.fa`
    
  + Quantification related inputs
    + Labelling method. Do Not use this option if you have label-free data, include the reference channel(s)
      for calculating relative peptide intensities. Possible options, tmtpro, tmt10plex, tmt6plex, tmt2plex, itraq8plex, itraq4plex.
      Multiple reference channels are averaged, multiple sets are separated by a space and must match the 
      `--mzmldef` parameter.
    
    `--isobaric 'set01:tmt10plex:130C:131 set02:tmt10plex:126'

  + Post-search processing inputs
  
    + map the genomic positions of VarDB peptides require the annotation GTF file of VarDB
    
    `--gtf /path/tovardb.gtf`

    + Canonical protein FASTA for catching canonical proteins and BLAST
    ```
    --blastdb /path/to/Uniprot.Ensembl.GENCODE.proteins.fa # Here we use latest uniprot, ensembl, gencode annotated protein sequences 
    --knownproteins /path/to/Homo_sapiens.GRCh38.pep.all.fa # a prefiltering known proteins DB needed to remove known peptides during class FDR calculation 
    ```
   
    + Annovar peptide annotation program location
    
    `--annovar_dir path/to/annovar # Downloaded before running, due to licensing`

    + Bigwig files for phastCons/PhyloCSF:
    `--bigwigs  /path/to/bigwigs/`

    + Mark novel peptides which can be explained by nsSNPs
    
    `--snpfa /path/to/MSCanProVar_ensemblV79.fa  # CanProVar annotated peptide sequences derived from known nsSNPs`   
    
    + Genome FASTA to BLAT against to find potential peptides mapped to multiple genomic locations
    
    `--genome /path/to/hg19.fa` # use hg19.fa.masked version if you don't want to consider repeated regions.
  
  + When using the VarDB database

    + SNP and COSMIC databases, required to map genomic positions of single amino acid variant peptides.
    ```
    --dbsnp /path/to/SNP142CodingDbSnp.txt # a text file containing genomic coordinates of coding SNPs
    --cosmic /path/to/CosmicMutantExport.tsv # a text file containing genomic coordinates of mutations
    ```
    
    + __Optional__: RNA-Seq BAM and BAI files (in same directory) for reads support in detected novel coding regions. 
    
    `--bamfiles '/path/to/*.bam'
  
  + Nextflow command option:
    + Use `-profile` option to run it in locally or submit it in slurm or sge system.
    ```
    -profile ## options are standard and testing. Options and cpus allocated can be re-defined in nextflow.config file.
    -resume  ## use it to resume the jobs from the last stopped process.
    ```
  + Nextflow configuration
    + Define CPU resources for specific processes in `configuration/base.config`
   

### Prepare once

  + Create account at [sanger](http://cancer.sanger.ac.uk/cosmic/help/download) for COSMIC database
  + [Register](http://annovar.openbioinformatics.org/en/latest) for download of annovar
  + Download SNP data from the [UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=661199271_5BEJQ6aAEOgRhkgNqBRFQQhTW05G&clade=mammal&org=&db=hg19&hgta_group=varRep&hgta_track=snp142Common&hgta_table=snp142CodingDbSnp&hgta_regionType=genome&position=&hgta_outputType=primaryTable&hgta_outFileName=snp142CodingDbSnp.txt)
  
```
# Get this repo
git clone https://github.com/lehtiolab/proteogenomics-analysis-workflow
cd proteogenomics-analysis-workflow

# Get Annovar
cd /path/to/your/annovar
wget __link_you_get_from_annovar__
tar xvfz annovar.latest.tar.gz
# This creates a folder with annotate_variation.pl and more files, to be passed to the pipeline with --annovar_dir


# Download bigwigs, this can take some time
cd /path/to/your/bigwigs  # this dir will be passed to the pipeline with --bigwigs
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons.bw 
wget https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg19/latest/PhyloCSF+0.bw
wget https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg19/latest/PhyloCSF+1.bw
wget https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg19/latest/PhyloCSF+2.bw
wget https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg19/latest/PhyloCSF-0.bw
wget https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg19/latest/PhyloCSF-1.bw
wget https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg19/latest/PhyloCSF-2.bw

# In the meantime, download and extract varDB data (Fasta, GTF, BlastP, SNP Fasta) to a good spot
wget -O varDB_data.tar.gz https://ndownloader.figshare.com/files/13358006 
tar xvfz varDB_data.tar.gz

# Get the hg19 masked genome sequence
wget hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFaMasked.tar.gz
tar xvfz chromFaMasked.gz
for chr in {1..22} X Y M; do cat chr$chr.fa.masked >> hg19.chr1-22.X.Y.M.fa.masked; done

# Download ENSEMBL database
wget ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
gunzip Homo_sapiens.GRCh38.pep.all.fa.gz

# Get the COSMIC database
sftp 'your_email_address@example.com'@sftp-cancer.sanger.ac.uk
# Download the data
sftp> get cosmic/grch37/cosmic/v81/CosmicMutantExport.tsv.gz
sftp> exit
# Extract COSMIC data
tar xvfz CosmicMutantExport.tsv.gz
```

### Analyse your mzML files with VarDB
Example command to search TMT 10-plex labelled data in docker
Remove  `--isobaric` parameter if you have label-free data.
```
nextflow run main.nf --tdb /path/to/VarDB.fasta \
  --mzmldef spectra_file_list.txt \
  --activation hcd \
  --isobaric 'set01:tmt10plex:131 set02:tmt10plex:131' 'set03:tmt10plex:127N' \
  --gtf /path/to/VarDB.gtf \
  --mods /path/to/tmt_mods.txt \
  --knownproteins /path/to/Homo_sapiens.GRCh38.pep.all.fa \
  --blastdb /path/to/UniProteome+Ensembl94+GENCODE24.proteins.fasta \
  --cosmic /path/to/CosmicMutantExport.tsv \
  --snpfa /path/to/MSCanProVar_ensemblV79.filtered.fasta \
  --genome /path/to/hg19.chr1-22.X.Y.M.fa \
  --dbsnp /path/to/snp142CodingDbSnp.txt \
  --annovar_dir /path/to/your/annovar \
  --bigwigs /path/to/your/bigwigs \
  --bamfiles /path/to/\*.bam \
  --outdir /path/to/results \
  -profile standard,docker # replace docker with singularity if needed
```
