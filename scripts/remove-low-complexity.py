# remove sequences with very low complexity
#
# needle-py needle/scripts/remove-low-complexity.py \
#   --forget-original \
#   `tangle-py tangle/scripts/defaults.py -m area_genomics_dir`/*


import shutil
import itertools
from pathlib import Path


def uniq_aa(v):
    lookup = {c:1 for c in v}
    return len(lookup.keys())


def filter_low_complexity(data_dict):
    return {k:v for k,v in data_dict.items() if uniq_aa(v) > 8}


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

        filtered = filter_low_complexity(protein_sequences)
        print(faa_fn, "filtered from", len(protein_sequences), "to", len(filtered))

        if not args.forget_original:
            orig_tsv_fn = tsv_fn+".orig"
            shutil.copy(tsv_fn, orig_tsv_fn)
            orig_faa_fn = faa_fn+".orig"
            shutil.copy(faa_fn, orig_faa_fn)

        write_fasta_from_dict(filtered, faa_fn)
        filtered_rows = [row for row in rows if row["target_accession"] in filtered]
        DetectedTable.write_tsv(tsv_fn, filtered_rows)
