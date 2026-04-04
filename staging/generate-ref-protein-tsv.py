import os
import csv
import argparse
import contextlib
from needle.duckdb import Unassigned
from needle.gff import parse_gff_to_hits
from scripts.defaults import DefaultPath
from needle.match import ProteinsTSV

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Generates protein TSV for classified NCBI RefSeq proteins.")
    ap.add_argument("module_id")
    ap.add_argument("output_protein_tsv", help="Output Protein TSV file")
    ap.add_argument("output_name_tsv", help="Output Protein TSV file")
    ap.add_argument("--overwrite", help="Overwrite old data", action="store_true", default=False)
    ap.add_argument("--genome-accession", help="Limit to these genome accession(s), can be a file", default=None)
    args = ap.parse_args()

    genome_accession_filter = []
    if args.genome_accession:
        if os.path.exists(args.genome_accession):
            with open(args.genome_accession, "r") as f:
                for acc in f.readlines():
                    acc = acc.strip()
                    genome_accession_filter.append(acc)
        else:
            genome_accession_filter.append(args.genome_accession)

    Unassigned.init(args.module_id)
    unassigned = Unassigned()
    protein_genome_accs = unassigned.protein_genome_accessions()
    genome_accessions = set([genome_acc for prot_acc, genome_acc in protein_genome_accs])
    protein_accessions = set([prot_acc for prot_acc, genome_acc in protein_genome_accs])

    genome_accessions = {acc for acc in genome_accessions \
      if os.path.exists(DefaultPath.ncbi_genome_gff(acc)) and (not genome_accession_filter or acc in genome_accession_filter)}

    print(f"{len(genome_accessions)} genome accessions")
    print(f"{len(protein_accessions)} protein accessions")
    name_tsv_fieldnames = ["protein_accession", "name"]

    if args.overwrite:
        print("overwriting previous data")
        with contextlib.suppress(FileNotFoundError):
            os.remove(args.output_protein_tsv)
            os.remove(args.output_name_tsv)

        with open(args.output_name_tsv, "w") as f:
            writer = csv.DictWriter(f, fieldnames=name_tsv_fieldnames, delimiter='\t')
            writer.writeheader()
    else:
        print("appending to previous data")

    for gac in genome_accessions:
        print(gac)
        gff_path = DefaultPath.ncbi_genome_gff(gac)
        proteins = parse_gff_to_hits(gff_path)
        proteins = [p for p in proteins if p.query_accession in protein_accessions]
        ProteinsTSV.append_to_tsv(args.output_protein_tsv, proteins, gac)

        with open(args.output_name_tsv, "a") as f:
            writer = csv.DictWriter(f, fieldnames=name_tsv_fieldnames, delimiter='\t')
            for p in proteins:
                if p._product_name is not None:
                    writer.writerow(dict(
                        protein_accession=p.query_accession,
                        name=p._product_name
                    ))
