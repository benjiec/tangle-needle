from tangle.sequence import read_fasta_as_dict
from needle.detect import get_aa_sequences
from Bio.Seq import Seq
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("fna_file")
ap.add_argument("contig")
ap.add_argument("start_1b", type=int)
ap.add_argument("stop_1b", type=int)
args = ap.parse_args()

d = read_fasta_as_dict(args.fna_file)
contig_sequence = d[args.contig]

if args.start_1b <= args.stop_1b:
    strand = 1
else:
    strand = -1

fragments = get_aa_sequences(
    args.contig,
    contig_sequence,
    min(args.start_1b, args.stop_1b),
    max(args.start_1b, args.stop_1b),
    strand,
    win=abs(args.start_1b-args.stop_1b)+1*2,
    win_overlap=0
)

if strand > 0:
    fragments = sorted(fragments, key=lambda t: t[1])
else:
    fragments = sorted(fragments, key=lambda t: t[2])[::-1]

for _,s,e,a in fragments:
    print(f">{args.contig}_{s}_{e}")
    print(a)
