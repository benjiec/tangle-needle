#!/usr/bin/env python3
import argparse

from needle.detect import Results
from needle.match import group_matches
from needle.clean import hmm_clean
from needle.expand import hmm_expand
from needle.hits import export_protein_hits
from needle.hmm import HMMCollection


def build_thresholds_dict(thresholds_file, threshold_frac):
    from heap.ko import KOThresholdTable
    from tangle.models import CSVSource

    src = CSVSource(KOThresholdTable, thresholds_file)
    values = src.values()

    thresholds = {}
    for row in values:
        t = row["threshold"]
        if t < 0:
           t = row["alen"]
        thresholds[row["model"]] = threshold_frac*t
    return thresholds


def main():
    parser = argparse.ArgumentParser(description="Export protein matches from detection results TSV.")
    parser.add_argument("hmm_file", help="HMM file to use for improve search results")
    parser.add_argument("fna_file")
    parser.add_argument("hmm_search_results_tsv")
    parser.add_argument("output_tsv_fn")
    parser.add_argument("output_faa_fn")
    parser.add_argument("--query-database-name", required=True)
    parser.add_argument("--append", action="store_true", default=False)
    parser.add_argument("--cpus", type=int, default=None)
    parser.add_argument("--threshold-file")
    parser.add_argument("--threshold-frac", default=0.8, type=float)
    args = parser.parse_args()

    res = Results(args.hmm_search_results_tsv, target_fasta_path=args.fna_file)
    protein_matches = group_matches(res.matches())
    accession_ids = [pm.query_accession for pm in protein_matches]
    hmm_collection = HMMCollection(args.hmm_file, accession_ids)

    try:
	# this step will further search around detected sequences, but now with
	# fewer queries (just 3 frame translations around the locus) and one
	# profile. because hmmsearch is used, the fewer sequences mean a weaker
	# match, perhaps in between previous matched fragments, may now match
	# with lower cond e-value (during detection e-value is used, not cond
	# e-value) and not be thrown away. this is helpful in filling in more
	# conserved or smaller fragments between and beyond already detected
	# fragments.

        thresholds = None

        if args.threshold_file:
            thresholds = build_thresholds_dict(args.threshold_file, args.threshold_frac)

        protein_matches = hmm_expand(
            protein_matches,
            res._target_sequences_by_accession,
            hmm_collection,
            thresholds = thresholds,
            cpus = args.cpus
        )

	# the cleaning step attempts to remove overlaps by finding the best
	# overlap (by how well the overlap matches to the HMM profile) between
	# every two fragment.

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
