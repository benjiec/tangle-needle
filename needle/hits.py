from tangle.detected import DetectedTable
from tangle import unique_batch


def write_fasta_record(f, pm) -> None:
    pid = pm.protein_hit_id
    seq = pm.collated_protein_sequence
    f.write(f">{pid}\n")
    f.write(seq + "\n")


def match_to_detected_row(protein_hit_id, match, genome_accession, batch):
    return dict(
        detection_type="model",  # target coordinates relative to target_model not target_accession
        detection_method="hmm",
        batch=batch,
        query_database=genome_accession,
        query_type="contig",
        query_accession=match.target_accession,
        target_database=genome_accession,
        target_type="protein",
        target_accession=protein_hit_id,
        target_model=match.query_accession,
        query_start=match.target_start,
        query_end=match.target_end,
        target_start=match.query_start,
        target_end=match.query_end
    )


def write_tsv(tsv_path, protein_hits, genome_accession, append):
    batch = unique_batch()
    rows = [match_to_detected_row(hit.protein_hit_id, m, genome_accession, batch) for hit in protein_hits for m in hit.matches]
    DetectedTable.write_tsv(tsv_path, rows, append = append)


def export_protein_hits(
    genome_accession,
    protein_hits,
    proteins_aa_path,
    proteins_tsv_path,
    append = False
) -> None:

    protein_hits = [pm for pm in protein_hits if pm.can_produce_single_sequence()]

    with open(proteins_aa_path, "at" if append else "wt") as f_prot:
        for pm in protein_hits:
            write_fasta_record(f_prot, pm)

    write_tsv(proteins_tsv_path, protein_hits, genome_accession, append)
