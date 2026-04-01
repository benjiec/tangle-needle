# Needle

Looking for a needle in a haystack...


## Setup

### HMM

Install hmmer3 package: e.g. on MacOS run `brew install hmmer`

Download the following profiles

  * KEGG KO profile HMMs: `https://www.genome.jp/ftp/db/kofam/profiles.tar.gz`
    * Then concatenate all the profiles together: `cat profiles/*.hmm > ko.hmm`
    * Run `hmmpress ko.hmm`
    * Also create `ko_thresholds.tsv` from KEGG FTP site `https://www.genome.jp/ftp/db/kofam/ko_list.gz`
      * Filter away rows without a threshold, using `grep -v "\-\t\-" ko_thresholds.txt`

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

Run these two commands to search for proteins matching HMM profile in a genomic FASTA file

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
