Integrated proteogenomics analysis workflow
==============

This is a workflow to identify, curate, and validate variant and novel peptides from MS proteomics spectra data, using the VarDB database. VarDB combines entries from COSMIC, PGOHUM, CanProVar and lncipedia. The pipeline takes mzML spectra files as input. The workflow is powered by [Nextflow](https://nextflow.io) and runs in [Docker](https://docker.com) containers.

Searches are run using [MSGF+](https://omics.pnl.gov/software/ms-gf) on 12 threads (adjust as you see fit) on a separate target and decoy databases which are then passed to [Percolator](http://percolator.ms) for statistical evaluation.


### Requirements

  + Linux system with
  + [Docker](https://docker.io)
  + [Nextflow](https://nextflow.io)
  + [Git](https://git-scm.com)


### Pipeline inputs

  + mzML files containing MS data
  
    `--mzmls /path/to/\*.mzML  # mind the backslash before *`
  + BAM and BAI files (in same directory) from RNASeq experiment
  
    `--bamfiles /path/to/\*.bam  # mind the backslash`
  + Target, decoy FASTA and GTF of VarDB

    `--tdb /path/to/vardb.fa --ddb /path/to/decoy_vardb.fa --gtf /path/tovardb.gtf`

  + Canonical protein FASTA for catching canonical proteins and BLAST
    `--knownproteins /path/to/Uniprot.Ensembl.RefSeq.GENCODE.proteins.fa`

  + SNP and COSMIC databases

    `--snpfa /path/to/SNPdb.fa` --dbsnp /path/to/SNP142CodingDbSnp.txt`
    `--cosmic /path/to/CosmicMutantExport.tsv`

  + Genome Masked FASTA to BLAT against

    `--genome /path/to/hg19.fa.masked`
  + Percolator pool size (fractionated data can be batched before percolator, this 
    specifies the amount of fractions in a batch. Batches will be merged and
    FDR re-evaluated after percolator

    `--ppoolsize 8  # default`
  + Isobaric quantification type used and activation (leave out for no isobaric quant)

    ```
    # default is NO isobaric quant. Use below or tmt6plex, tmt2plex, itraq8plex, itraq4plex 
    --isobaric tmt10plex
    --activation hcd  # default else use cid, etd
    ```

### Prepare once

  + Create account at [sanger](http://cancer.sanger.ac.uk/cosmic/help/download) for COSMIC database
  + [Register](http://annovar.openbioinformatics.org/en/latest) for download of annovar
  + Download SNP data from the [UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=654845801_AIfwaTHVOpBosVlaTdk1QGgcQYrZ&clade=mammal&org=Human&db=hg38&hgta_group=varRep&hgta_track=snp142Common&hgta_table=snp142CodingDbSnp&hgta_regionType=genome&position=chr1%3A11102837-11267747&hgta_outputType=primaryTable&hgta_outFileName=snp142CodingDbSnp.txt)
  
```
# Get this repo
git clone https://github.com/proteogenomics-analysis-workflow
cd proteogenomics-analysis-workflow
# Get Annovar
wget __link_you_get_from_annovar__
tar xvfz annovar.latest.tar.gz

# Create containers
docker build -f spectrumAI_Dockerfile -t spectrumai .
docker build -f annovar_Dockerfile -t annovar .
docker build -f pgpython_Dockerfile -t pgpython .  # downloads bigwig files, takes a long time

# In the meantime, download and extract varDB data (Fasta, GTF, BlastP, SNP Fasta to a good spot
wget http://lehtiolab.se/Supplementary_Files/varDB_data.tar.gz
tar xvfz varDB_data.tar.gz

# Get the hg19 masked genome sequence
wget hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFaMasked.tar.gz
tar xvfz chromFaMasked.gz

# Get the COSMIC database
sftp 'your_email_address@example.com'@sftp-cancer.sanger.ac.uk
# Download the data (NB you may want to download other versions)
sftp> get cosmic/grch37/cosmic/v75/CosmicMutantExport.tsv.gz
sftp> exit
# Extract COSMIC data
tar xvfz CosmicMutantExport.tsv.gz
```

### Analyse your mzML files

```
nextflow run ipaw.nf --tdb /path/to/VarDB.fasta --ddb /path/to/decoy_VarDB.fasta \ 
  --mzmls /path/to/\*.mzML --gtf /path/to/VarDB.gtf \
  --knownproteins /path/to/UniProteome+Ensembl87+refseq+GENCODE24.proteins.fasta \
  --snpfa /path/to/MSCanProVar_ensemblV79.filtered.fasta \
  --genome /path/to/hg19.chr1-22.X.Y.M.fa.masked \
  --dbsnp /path/to/snp142CodingDbSnp.txt
  --bamfiles /path/to/\*.bam --isobaric tmt10plex
```
