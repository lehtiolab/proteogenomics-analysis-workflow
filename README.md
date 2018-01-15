Proteogenomics analysis workflow
==============

This is a workflow to identify, curate, and validate variant and novel peptides from MS proteomics spectra data. It is powered by [Nextflow](https://nextflow.io) and runs in [Docker](https://docker.com) containers.

To run:

 + clone this repo:
```
git clone https://github.com/proteogenomics-analysis-workflow
```

 + install [Nextflow](https://nextflow.io)
`sudo apt install nextflow  # this assumes an Ubuntu distro`

 + Run and sit back
`nextflow run paw.nf --mzmls /path/to/mzmlfiles --tdb /path/to/target_db --ddb /path/to/decoydb --gtf /path/to/gtf --bed /path/to/bedfiles`

This assumes FASTA data is from a variant database called VarDB, which contains entries from COSMIC, PGOHUM, CanProVar and lncipedia.
