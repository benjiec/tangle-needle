# Needle

Looking for a needle in a haystack...


## Setup

### HMM

Install hmmer3 package: e.g. on MacOS run `brew install hmmer`

Download the following profiles

  * KEGG KO profile HMMs: `https://www.genome.jp/ftp/db/kofam/profiles.tar.gz`
    * Then concatenate all the profiles together: `cat profiles/*.hmm > ko.hmm`
      * Remove all entries that are RNA or small RNAs: make a list of KOs to include, then use `hmmfetch -f` to create a new `.hmm` file
      * Run `hmmfetch --index ko.hmm`

Put these files in the same directory then set the following environment variable

```
HMM_DB_DIR=</host/hmm_dir>
```

## Hosting Docker Images on Google Cloud

This repo

```
docker build --platform linux/amd64 -t us-east1-docker.pkg.dev/needle-489321/tangle-docker/needle:latest .
docker push us-east1-docker.pkg.dev/needle-489321/tangle-docker/needle:latest
```

And, to use these images on Google Cloud, make sure everything under HMM_DB_DIR
are synced to Google Cloud storage, into a bucket.


## Protein Detection from HMM Profiles

Run these two commands to search for proteins matching HMM profile in a genomic
FASTA file. Note, the second command speeds up significantly with `.ssi` file
generated from `hmmfetch --index`.

```
python3 scripts/hmmsearch-genome.py \
  <hmm-file> <genomic-fasta-file> hmm-search-output.tsv
python3 scripts/export-protein-results.py \
  --query-database-name <genome-accession> \
  <hmm-file> <genomic-fasta-file> hmm-search-output.tsv \
  output_proteins.tsv output_proteins.faa
```

The final outputs are two files: a TSV file that enumerates the fragments on
the contigs (this is like a GFF, in the tangle DetectedTable format), and a
protein FASTA file.

To prepare and submit this job to run on Google Cloud, use the following script
to create a run directory under `runs` (or whatever value for `--run-dir`), and
then follow instructions in the README file in that run directory.

```
python3 gcloud/hmm-detect/setup.py \
  --genome-accession GCF_002042975.1 \
  --run-dir=runs \
  ncbi-downloads/ncbi_dataset/data/GCF_002042975.1/GCF_002042975.1_ofav_dov_v1_genomic.fna
```
