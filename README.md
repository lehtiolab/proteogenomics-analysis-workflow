Integrated proteogenomics analysis workflow
==============

This is a workflow to identify, curate, and validate variant and novel peptides from MS proteomics spectra data, using the VarDB database. VarDB combines entries from COSMIC, PGOHUM, CanProVar and lncipedia. The pipeline takes mzML spectra files as input. The workflow is powered by [Nextflow](https://nextflow.io) and runs in [Docker](https://docker.com) containers.

Searches are run using [MSGF+](https://omics.pnl.gov/software/ms-gf) on 12 threads (adjust as you see fit) on a concatenated target and decoy databases which are then passed to [Percolator](http://percolator.ms) for statistical evaluation.

Please cite this paper when you have used the workflow for publications.

Zhu Y, Orre LM, Johansson HJ, Huss M, Boekel J, Vesterlund M, Fernandez-Woodbridge A, Branca RMM, Lehtio J: Discovery of coding regions in the human genome by integrated proteogenomics analysis workflow. Nat Commun 2018, 9(1):903.  [PMID: 29500430](https://www.ncbi.nlm.nih.gov/pubmed/29500430)

![workflow image](https://github.com/lehtiolab/proteogenomics-analysis-workflow/blob/master/images/workflow.png)

### Requirements

  + Linux system with
  + [Docker](https://docker.io)
  + [Nextflow](https://nextflow.io)
  + [Git](https://git-scm.com)


### Pipeline inputs

  + mzML files containing MS data
  
    `--mzmls /path/to/\*.mzML  # mind the backslash before *`
  + __Optional__: BAM and BAI files (in same directory) from RNASeq experiment
  
    `--bamfiles /path/to/\*.bam  # mind the backslash`
  + FASTA and GTF of VarDB

    `--tdb /path/to/vardb.fa --gtf /path/tovardb.gtf`

  + Modification file for MSGF+. Default file is for TMT samples. [Here is an example.](https://bix-lab.ucsd.edu/download/attachments/13533355/Mods.txt?version=2&modificationDate=1358975546000)
    `--mods Mods.txt`

  + Canonical protein FASTA for catching canonical proteins and BLAST
    `--blastdb /path/to/Uniprot.Ensembl.RefSeq.GENCODE.proteins.fa`
    `--knownproteins /path/to/Homo_sapiens.GRCh38.pep.all.fa`

  + SNP and COSMIC databases

    `--snpfa /path/to/SNPdb.fa --dbsnp /path/to/SNP142CodingDbSnp.txt`
    `--cosmic /path/to/CosmicMutantExport.tsv`

  + Genome Masked FASTA to BLAT against

    `--genome /path/to/hg19.fa.masked`
  + Isobaric quantification type used and activation (leave out for no isobaric quant)

    ```
    # default is NO isobaric quant (label-free). Don't use this option unless for tmt10plex, tmt6plex, tmt2plex, itraq8plex, itraq4plex
    --isobaric tmt10plex
    --activation hcd  # default else use cid, etd
    ```

### Prepare once

  + Create account at [sanger](http://cancer.sanger.ac.uk/cosmic/help/download) for COSMIC database
  + [Register](http://annovar.openbioinformatics.org/en/latest) for download of annovar
  + Download SNP data from the [UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=661199271_5BEJQ6aAEOgRhkgNqBRFQQhTW05G&clade=mammal&org=&db=hg19&hgta_group=varRep&hgta_track=snp142Common&hgta_table=snp142CodingDbSnp&hgta_regionType=genome&position=&hgta_outputType=primaryTable&hgta_outFileName=snp142CodingDbSnp.txt)
  
```
# Get this repo
git clone https://github.com/proteogenomics-analysis-workflow
cd proteogenomics-analysis-workflow

# Get Annovar
cd dockerfiles
wget __link_you_get_from_annovar__
tar xvfz annovar.latest.tar.gz

# Download bigwigs
docker build -f pgpython_bigwigs -t pgpython_bigwigs .  # downloads bigwig files, takes a long time

# Create pipeline containers
docker build -f spectrumAI -t spectrumai .
docker build -f annovar_Dockerfile -t annovar .
docker build -f pgpython -t pgpython . 
cd ..

# In the meantime, download and extract varDB data (Fasta, GTF, BlastP, SNP Fasta) to a good spot
wget http://lehtiolab.se/Supplementary_Files/varDB_data.tar.gz
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
# Download the data (NB version 71 currently works with the mapping script)
sftp> get cosmic/grch37/cosmic/v71/CosmicMutantExport.tsv.gz
sftp> exit
# Extract COSMIC data
tar xvfz CosmicMutantExport.tsv.gz
```

### Analyse your mzML files
Example command to search TMT 10-plex labelled data.
Remove  `--isobaric tmt10plex`  if you have label-free data.
```
nextflow run ipaw.nf --tdb /path/to/VarDB.fasta \ 
  --mzmls /path/to/\*.mzML --gtf /path/to/VarDB.gtf \
  --mods /path/to/mods.txt \
  --knownproteins /path/to/Homo_sapiens.GRCh38.pep.all.fa \
  --blastdb /path/to/UniProteome+Ensembl87+refseq+GENCODE24.proteins.fasta \
  --cosmic /path/to/CosmicMutantExport.tsv \
  --snpfa /path/to/MSCanProVar_ensemblV79.filtered.fasta \
  --genome /path/to/hg19.chr1-22.X.Y.M.fa.masked \
  --dbsnp /path/to/snp142CodingDbSnp.txt \
  --bamfiles /path/to/\*.bam --isobaric tmt10plex \
  --outdir /path/to/results
```
