from tangle.sequence import read_fasta_as_dict
from Bio.Seq import Seq
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("fna_file")
ap.add_argument("contig")
ap.add_argument("start_1b", type=int)
ap.add_argument("stop_1b", type=int)
ap.add_argument("name")
args = ap.parse_args()

d = read_fasta_as_dict(args.fna_file)
contig_sequence = d[args.contig]
sequence = contig_sequence[args.start_1b-1:args.stop_1b]

print(f">{args.name}")
print(sequence)
