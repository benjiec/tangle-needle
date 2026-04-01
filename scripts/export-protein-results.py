#!/usr/bin/env python3
import argparse

from needle.detect import Results
from needle.match import group_matches, export_protein_hits
from needle.hits import hmm_expand, hmm_clean
from needle.hmm import HMMCollection
from scripts.defaults import DefaultPath

def main():
    parser = argparse.ArgumentParser(description="Export protein matches from detection results TSV.")
    parser.add_argument("results_tsv", help="Detection results TSV (Results.PRODUCER headers)")
    parser.add_argument("genome_accession", help="Genome accession")
    parser.add_argument("hmm_file", help="HMM file to use for improve search results")
    parser.add_argument("output_dir", help="Path of output directory")
    parser.add_argument("--fna", help="HMM accession to export", default=None)
    parser.add_argument("--hmm", help="HMM accession to export", default=None)
    parser.add_argument("--no-output", help="Skip output", action='store_true', default=False)
    args = parser.parse_args()

    if args.fna:
        target_fasta = args.fna
    else:
        target_fasta = DefaultPath.ncbi_genome_fna(args.genome_accession)

    res = Results(args.results_tsv, target_fasta_path=target_fasta)
    protein_matches = group_matches(res.matches())
    accession_ids = [pm.query_accession for pm in protein_matches]
    if args.hmm:
        protein_matches = [m for m in protein_matches if m.query_accession == args.hmm]
        accession_ids = [args.hmm]
    hmm_collection = HMMCollection(args.hmm_file, accession_ids)

    try:

        # use HMM to refine protein search at detected locus
	# this step helps because we are using conditinoal e-value rather than
	# independent e-value, for matches, which would be more sensitive given
	# we already decided the locus is where protein is
        protein_matches = hmm_expand(protein_matches, res._target_sequences_by_accession, hmm_collection)

        pre_filter = len(protein_matches)
        protein_matches = [m for m in protein_matches if m.can_collate()]
        print("filter by collatable", pre_filter, "=>", len(protein_matches))
        cleaned_protein_matches = hmm_clean(protein_matches, hmm_collection)

        if args.no_output is False:
            export_protein_hits(
                args.genome_accession,
                cleaned_protein_matches,
                args.output_dir+"/protein_detected.faa",
                args.output_dir+"/protein_detected.tsv"
            )

    finally:
        hmm_collection.clean()


if __name__ == "__main__":
    main()
