#!/usr/bin/env python3
import argparse

from needle.detect import Results
from needle.match import group_matches
from needle.hits import hmm_expand, hmm_clean, export_protein_hits
from needle.hmm import HMMCollection

def main():
    parser = argparse.ArgumentParser(description="Export protein matches from detection results TSV.")
    parser.add_argument("hmm_file", help="HMM file to use for improve search results")
    parser.add_argument("fna_file")
    parser.add_argument("hmm_search_results_tsv")
    parser.add_argument("output_tsv_fn")
    parser.add_argument("output_faa_fn")
    parser.add_argument("--query-database-name", required=True)
    parser.add_argument("--append", action="store_true", default=False)
    args = parser.parse_args()

    res = Results(args.hmm_search_results_tsv, target_fasta_path=args.fna_file)
    protein_matches = group_matches(res.matches())
    accession_ids = [pm.query_accession for pm in protein_matches]
    hmm_collection = HMMCollection(args.hmm_file, accession_ids)

    try:
	# skipping this now: seems like normal detection through the contig
	# works ok and the expansion logic here is expensive w/o getting finer
	# results.
        #
        # protein_matches = hmm_expand(protein_matches, res._target_sequences_by_accession, hmm_collection)

        pre_filter = len(protein_matches)
        protein_matches = [m for m in protein_matches if m.can_collate()]
        print("filter by collatable", pre_filter, "=>", len(protein_matches))
        # cleanup overlaps between matches
        cleaned_protein_matches = hmm_clean(protein_matches, hmm_collection)

        export_protein_hits(
            args.query_database_name,
            cleaned_protein_matches,
            args.output_faa_fn,
            args.output_tsv_fn,
            append=args.append
        )

    finally:
        hmm_collection.clean()


if __name__ == "__main__":
    main()
