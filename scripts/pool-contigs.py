# pool multiple genome fasta files together, adding the genome accession to each contig name

import argparse
from pathlib import Path
from tangle.sequence import read_fasta_as_dict, write_fasta_from_dict

ap = argparse.ArgumentParser()
ap.add_argument("pooled_fna")
ap.add_argument("genome_fnas", nargs="+")
args = ap.parse_args()

for genome_fna in args.genome_fnas:
    contigs = read_fasta_as_dict(genome_fna)
    genome_accession = Path(genome_fna).parent.name
    print(f"Using {genome_accession} as accession for {genome_fna}")
    converted = {}
    for acc, seq in contigs.items():
        converted[f"{genome_accession}|{acc}"] = seq
    write_fasta_from_dict(converted, args.pooled_fna, append=True)
