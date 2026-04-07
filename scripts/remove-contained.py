# remove duplicated entries based on sequence

import itertools
import ahocorasick

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
    ap.add_argument("tsv_fn")
    ap.add_argument("fasta_fn")
    ap.add_argument("filtered_tsv_fn")
    ap.add_argument("filtered_fasta_fn")

    args = ap.parse_args()

    source = CSVSource(DetectedTable, args.tsv_fn)
    rows = source.values()
    protein_sequences = read_fasta_as_dict(args.fasta_fn)

    filtered = filter_contained_sequences(rows, protein_sequences)
    print("filtered from", len(protein_sequences), "to", len(filtered))

    write_fasta_from_dict(filtered, args.filtered_fasta_fn)
    filtered_rows = [row for row in rows if row["target_accession"] in filtered]
    DetectedTable.write_tsv(args.filtered_tsv_fn, filtered_rows)
