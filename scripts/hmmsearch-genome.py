import argparse
from tangle.sequence import read_fasta_as_dict
from needle.detect import hmm_search_genome, Results

ap = argparse.ArgumentParser()
ap.add_argument("hmm_file")
ap.add_argument("fna_file")
ap.add_argument("output_file")
ap.add_argument("--append", action="store_true", default=False)
ap.add_argument("--cpus", type=int, default=None)
ap.add_argument("--target-accession", type=str, default=None)
ap.add_argument("--target-left", type=int, default=None)
ap.add_argument("--target-right", type=int, default=None)
args = ap.parse_args()

fna_file = args.fna_file
genomic_fasta = read_fasta_as_dict(fna_file)

hmm_rows = hmm_search_genome(
    args.hmm_file, genomic_fasta,
    target_accession = args.target_accession,
    target_left = args.target_left,
    target_right = args.target_right,
    cpus = args.cpus
)

detected = []
for row in hmm_rows:
    out = (
      row["query_accession"],
      row["target_accession"],
      row["dom_evalue"],
      '',
      row["hmm_from"],
      row["hmm_to"],
      row["ali_from"],
      row["ali_to"],
      row["matched_sequence"]
    )
    detected.append(out)

with open(args.output_file, "a" if args.append else "w") as f:
    if args.append is False:
        f.write("\t".join(Results.PRODUCER_HEADER)+"\n")
    for d in detected:
        assert len(d) == len(Results.PRODUCER_HEADER)
        f.write("\t".join([str(x) for x in d])+"\n")
