# demultiplex multiple genome's faa and TSV outputs, from using genome_accession|contig in fasta, to listing genome_accesss as query and target databases

import argparse
import duckdb
from tangle.sequence import read_fasta_as_dict, write_fasta_from_dict
from tangle.models import Schema, CSVSource
from tangle.detected import DetectedTable

ap = argparse.ArgumentParser()
ap.add_argument("tsv_fn")
ap.add_argument("faa_fn")
ap.add_argument("out_tsv_fn")
ap.add_argument("out_faa_fn")
args = ap.parse_args()

source = CSVSource(DetectedTable, args.tsv_fn)
schema = Schema("schema")
schema.add_table(source)
schema.duckdb_load()
db = duckdb.connect(':default:')
query = f"SELECT * FROM {schema.name}.{DetectedTable.name}"
rows = db.execute(query).fetchdf().to_dict('records')

proteins = read_fasta_as_dict(args.faa_fn)


def convert_accession(a):
    if "|" in a:
        splitted = a.split("|")
        return splitted[0], "|".join(splitted[1:])
    return None, a


for row in rows:
    query_db, query_acc = convert_accession(row["query_accession"])
    if query_db:
        row["query_database"] = query_db
        row["query_accession"] = query_acc

    target_db, target_acc = convert_accession(row["target_accession"])
    if target_db:
        row["target_database"] = target_db
        row["target_accession"] = target_acc

corrected_proteins = {}
for acc, sequence in proteins.items():
    _, new_acc = convert_accession(acc)
    corrected_proteins[new_acc] = sequence

write_fasta_from_dict(corrected_proteins, args.out_faa_fn)
DetectedTable.write_tsv(args.out_tsv_fn, rows)
