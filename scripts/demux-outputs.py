# demultiplex multiple genome's faa and TSV outputs, from using genome_accession|contig in fasta, to listing genome_accesss as query and target databases

import argparse
from pathlib import Path
from tangle.sequence import read_fasta_as_dict, write_fasta_from_dict
from tangle.models import CSVSource
from tangle.detected import DetectedTable

ap = argparse.ArgumentParser()
ap.add_argument("tsv_fn")
ap.add_argument("faa_fn")
ap.add_argument("demuxed_tsv_fn")
ap.add_argument("output_faa_dir")
ap.add_argument("--output_faa_file_name", default="proteins.faa")

args = ap.parse_args()

source = CSVSource(DetectedTable, args.tsv_fn)
rows = source.values()

proteins = read_fasta_as_dict(args.faa_fn)
proteins_by_db = {}

def convert_accession(a):
    if "|" in a:
        splitted = a.split("|")
        return splitted[0], "|".join(splitted[1:])
    return None, a

for row in rows:
    if row["target_accession"] not in proteins:
        print(f"warning: cannot find {row['target_accession']} in {args.faa_fn}")
        sequence = None
    else:
        sequence = proteins[row["target_accession"]]

    query_db, query_acc = convert_accession(row["query_accession"])
    if query_db:
        row["query_database"] = query_db
        row["query_accession"] = query_acc

    target_db, target_acc = convert_accession(row["target_accession"])
    if target_db:
        row["target_database"] = target_db
        row["target_accession"] = target_acc
        if target_db not in proteins_by_db:
            proteins_by_db[target_db] = {}
        if sequence:
            proteins_by_db[target_db][target_acc] = sequence
    else:
        print(f"warning: cannot determine a target_database value for {row['target_accession']}, skip writing to a fasta file")

DetectedTable.write_tsv(args.demuxed_tsv_fn, rows)

for db, seq_dict in proteins_by_db.items():
    parent_parent_dir = Path(args.output_faa_dir)
    parent_parent_dir.mkdir(exist_ok=True)
    parent_dir = parent_parent_dir / db
    parent_dir.mkdir(exist_ok=True)
    fasta_fn = parent_dir / args.output_faa_file_name
    write_fasta_from_dict(seq_dict, str(fasta_fn))
