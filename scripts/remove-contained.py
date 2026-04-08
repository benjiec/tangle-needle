# remove duplicated entries based on sequence

import shutil
import itertools
import ahocorasick
from pathlib import Path


def filter_contained_sequences(rows, data_dict):
    # deduplicate exact strings first to reduce trie size
    unique_to_keys = {}
    for k, v in data_dict.items():
        if v not in unique_to_keys:
            unique_to_keys[v] = k

    sorted_items = sorted(
      unique_to_keys.items(),
      key=lambda x: len(x[0]), reverse=True
    )

    # build the automaton with all unique sequences
    automaton = ahocorasick.Automaton()
    for seq, key in sorted_items:
        automaton.add_word(seq, (key, seq))
    automaton.make_automaton()
    print("done making automaton")

    # a locus is (query_accession, query_left, query_right, strand), where
    # query_accession is essentially the contig accession

    target_accession_locus = {}
    group_fn = lambda row: row["target_accession"]
    rows = sorted(rows, key=group_fn)
    for tacc, group in itertools.groupby(rows, group_fn):
        group = list(group)
        query_left = min([min(row["query_start"], row["query_end"]) for row in group])
        query_right = max([max(row["query_start"], row["query_end"]) for row in group])
        strand = 1 if group[0]["query_start"] < group[0]["query_end"] else -1
        locus = (group[0]["query_accession"], query_left, query_right, strand)
        target_accession_locus[group[0]["target_accession"]] = locus

    def in_locus(big, small):
        return \
            big[0] == small[0] and \
            big[3] == small[3] and \
            big[1] <= small[1] and \
            big[2] >= small[2]

    to_remove = set()
    for seq, key in sorted_items:
        if key in to_remove:
            continue
        seq_locus = target_accession_locus[key]

        # search 'seq' against the automaton to find any internal patterns
        for end_index, (found_contained_key, found_contained_seq) in automaton.iter(seq):
            # if the found pattern is not the string itself, it's a contained sequence
            if found_contained_key != key:
		# only remove that contained sequence if it's at the same
		# locus, this way we keep duplicated sequences from other parts
		# of the genome around
                contained_locus = target_accession_locus[found_contained_key]
                if in_locus(seq_locus, contained_locus):
                    to_remove.add(found_contained_key)

    # 4. Return dictionary excluding the contained keys
    return {k: v for k, v in data_dict.items() if k in unique_to_keys.values() and k not in to_remove}


if __name__ == "__main__":
    import argparse
    from tangle.sequence import read_fasta_as_dict, write_fasta_from_dict
    from tangle.models import CSVSource
    from tangle.detected import DetectedTable

    ap = argparse.ArgumentParser()
    ap.add_argument("genome_dirs", nargs="+")
    ap.add_argument("--forget-original", action="store_true", default=False)
    ap.add_argument("--demuxed-fasta-filename", default="proteins.faa")
    ap.add_argument("--demuxed-tsv-filename", default="proteins.tsv")
    args = ap.parse_args()

    for genome_dir in args.genome_dirs:
        genome_dir_path = Path(genome_dir)
        tsv_fn = str(genome_dir_path / args.demuxed_tsv_filename)
        faa_fn = str(genome_dir_path / args.demuxed_fasta_filename)

        source = CSVSource(DetectedTable, tsv_fn)
        rows = source.values()
        protein_sequences = read_fasta_as_dict(faa_fn)

        filtered = filter_contained_sequences(rows, protein_sequences)
        print(faa_fn, "filtered from", len(protein_sequences), "to", len(filtered))

        if not args.forget_original:
            orig_tsv_fn = tsv_fn+".orig"
            shutil.copy(tsv_fn, orig_tsv_fn)
            orig_faa_fn = faa_fn+".orig"
            shutil.copy(faa_fn, orig_faa_fn)

        write_fasta_from_dict(filtered, faa_fn)
        filtered_rows = [row for row in rows if row["target_accession"] in filtered]
        DetectedTable.write_tsv(tsv_fn, filtered_rows)
