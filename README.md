Proteogenomics analysis workflow
==============

This is a workflow to identify, curate, and validate variant and novel peptides from MS proteomics spectra data, using the VarDB database. VarDB combines entries from COSMIC, PGOHUM, CanProVar and lncipedia. The workflow is powered by [Nextflow](https://nextflow.io) and runs in [Docker](https://docker.com) containers.


### Requirements

  + Linux system with
  + Docker 
  + Nextflow
  + Git

### Input data

  + mzML files containing MS data `--mzmls /path/to/\*.mzML  # mind the backslash before *`
  + Target, decoy FASTA formatted VarDB files `--tdb /path/to/vardb.fa --ddb /path/to/decoy_vardb.fa`
  + GTF file of VarDB `--gtf /path/tovardb.gtf`
  + BlastP FASTA database `--blastdb  /path/to/Uniprot.Ensembl.RefSeq.GENCODE.proteins.fa`
  + SNP FASTA database `--snpdb /path/to/SNPdb.fa`
  + Genome Masked FASTA `--genome /path/to/hg19.fa.masked`
  + BAM and BAI files (in same directory) from RNASeq experiment `--bamfiles /path/to/\*.bam  # backslash to escape *`

### To run

  + Prepare once:
```
# Get this repo
git clone https://github.com/proteogenomics-analysis-workflow
cd proteogenomics-analysis-workflow

# Create containers
docker build -f pgpython_Dockerfile -t pgpython .  # downloads bigwig files, takes a long time
docker build -f spectrumAI_Dockerfile -t spectrumai .

# Download FASTA databases

```
  + Run pipeline:
```
sudo nextflow run paw.nf --tdb /path/to/VarDB.3frame.fasta --ddb /path/to/decoy_VarDB.3frame.fasta --mzmls /path/to/\*.mzML --gtf /path/to/VarDB.gtf --blastdb /path/to/UniProteome+Ensembl87+refseq+GENCODE24.proteins.fasta  --snpdb /path/to/MSCanProVar_ensemblV79.filtered.fasta --genome /path/to/hg19.chr1-22.X.Y.M.fa.masked  --bamfiles /path/to/\*.bam
```
