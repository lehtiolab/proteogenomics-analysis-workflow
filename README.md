Integrated proteogenomics analysis workflow
==============

This is a workflow to identify, curate, and validate variant and novel peptides from MS proteomics spectra data, using the VarDB database. VarDB combines entries from COSMIC, PGOHUM, CanProVar and lncipedia. The pipeline takes mzML spectra files as input. The workflow is powered by [Nextflow](https://nextflow.io) and runs in [Docker](https://docker.com) containers.


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

    `--snpfa /path/to/SNPdb.fa` --dbsnp /path/to/SNP142CodingDbSnp.txt
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

### To run

  + Prepare once:
  
```
# Get this repo
git clone https://github.com/proteogenomics-analysis-workflow
cd proteogenomics-analysis-workflow

# Create containers
docker build -f spectrumAI_Dockerfile -t spectrumai .
docker build -f pgpython_Dockerfile -t pgpython .  # downloads bigwig files, takes a long time

# In the meantime, download and extract varDB data (Fasta, GTF, BlastP, SNP Fasta
wget https://lehtiolab.se/files/vardb.tar.gz
tar xvfz vardb.tar.gz

# Download COSMIC, SNP DB
wget https:// snp, hg19
```

  + Analyse your mzML files

```
sudo nextflow run ipaw.nf --tdb /path/to/VarDB.fasta --ddb /path/to/decoy_VarDB.fasta \ 
  --mzmls /path/to/\*.mzML --gtf /path/to/VarDB.gtf \
  --knownproteins /path/to/UniProteome+Ensembl87+refseq+GENCODE24.proteins.fasta \
  --snpfa /path/to/MSCanProVar_ensemblV79.filtered.fasta --genome /path/to/hg19.chr1-22.X.Y.M.fa.masked \
  --bamfiles /path/to/\*.bam --isobaric tmt10plex
```
