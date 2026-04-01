import csv
from typing import Dict, List, Optional
from Bio.Seq import Seq

from tangle.sequence import read_fasta_as_dict
from .match import Match
from .sequence import extract_subsequence, extract_subsequence_strand_sensitive
from .sequence import to_dna_coordinate, compute_three_frame_translations
from .hmm import hmmsearch_sequence_dict


DOM_EVALUE_LIMIT = 0.01


class Results:
    # NCBI-style canonical headers (preferred)
    H_QSEQID = "qseqid"
    H_SSEQID = "sseqid"
    H_EVALUE = "evalue"
    H_PIDENT = "pident"
    H_QSTART = "qstart"
    H_QEND = "qend"
    H_SSTART = "sstart"
    H_SEND = "send"
    H_SSEQ = "sseq"

    # Producer header order (used by blast-genome.py). Keep minimal set we rely on.
    PRODUCER_HEADER = [
        H_QSEQID,
        H_SSEQID,
        H_EVALUE,
        H_PIDENT,
        H_QSTART,
        H_QEND,
        H_SSTART,
        H_SEND,
        H_SSEQ,
    ]

    # Raw BLAST outfmt fields (as configured in run_blast_search)
    RAW_OUTFMT_FIELDS = [
        "qseqid",
        "sseqid",
        "evalue",
        "pident",
        "qstart",
        "qend",
        "stitle",  # not included in output header, but present in raw
        "sseq",
        "sstart",
        "send",
    ]

    # Mapping from raw BLAST field names to our header names
    RAW_TO_HEADER = {
        "qseqid": H_QSEQID,
        "sseqid": H_SSEQID,
        "evalue": H_EVALUE,
        "pident": H_PIDENT,
        "qstart": H_QSTART,
        "qend": H_QEND,
        "sseq": H_SSEQ,
        "sstart": H_SSTART,
        "send": H_SEND,
    }
    # Inverse mapping: header -> raw key
    HEADER_TO_RAW = {v: k for k, v in RAW_TO_HEADER.items()}

    def __init__(
        self,
        results_tsv_path: str,
        target_fasta_path: Optional[str] = None
    ) -> None:
        self._results_tsv_path = results_tsv_path
        self._target_fasta_path = target_fasta_path
        self._matches: Optional[List[Match]] = None
        self._target_sequences_by_accession: Optional[Dict[str, str]] = None

    def matches(self) -> List[Match]:
        if self._matches is None:
            self._parse_once()
        # mypy: _matches is now set
        return self._matches or []

    # ----- Internal helpers -----

    def _parse_once(self) -> None:
        if self._target_fasta_path:
            self._target_sequences_by_accession = read_fasta_as_dict(self._target_fasta_path)

        with open(self._results_tsv_path, "r") as tsv_file:
            reader = csv.reader(tsv_file, delimiter="\t")
            header_row = next(reader, None)
            if header_row is None:
                self._matches = []
                return
            header_index = self._build_header_index(header_row)

            matches: List[Match] = []
            for row in reader:
                if not row or all(not cell for cell in row):
                    continue
                match = self._row_to_match(row, header_index)
                # Attach sequences if requested and available
                if match is not None:
                    if self._target_sequences_by_accession is not None:
                        seq = extract_subsequence_strand_sensitive(
                            self._target_sequences_by_accession.get(match.target_accession, None),
                            match.target_start,
                            match.target_end,
                        )
                        match.target_sequence = seq
                    matches.append(match)

            self._matches = matches

    def _build_header_index(self, header_row: List[str]) -> Dict[str, int]:
        # Normalize header names by stripping whitespace
        normalized = [h.strip() for h in header_row]
        index: Dict[str, int] = {name: i for i, name in enumerate(normalized)}

        # Build a mapping for canonical names only; sseq is optional
        mapping: Dict[str, Optional[int]] = {
            self.H_QSEQID: index.get(self.H_QSEQID),
            self.H_SSEQID: index.get(self.H_SSEQID),
            self.H_EVALUE: index.get(self.H_EVALUE),
            self.H_PIDENT: index.get(self.H_PIDENT),
            self.H_QSTART: index.get(self.H_QSTART),
            self.H_QEND: index.get(self.H_QEND),
            self.H_SSTART: index.get(self.H_SSTART),
            self.H_SEND: index.get(self.H_SEND),
            self.H_SSEQ: index.get(self.H_SSEQ),
        }

        # Ensure required fields exist
        required = [
            self.H_QSEQID,
            self.H_SSEQID,
            self.H_EVALUE,
            self.H_PIDENT,
            self.H_QSTART,
            self.H_QEND,
            self.H_SSTART,
            self.H_SEND,
        ]
        missing = [key for key in required if mapping.get(key) is None]
        if missing:
            raise ValueError(f"Missing required columns in results TSV: {', '.join(missing)}")

        # Convert Optional[int] to int where present
        return {k: v for k, v in mapping.items() if v is not None}

    def _row_to_match(self, row: List[str], header_index: Dict[str, int]) -> Optional[Match]:
        qacc = row[header_index[self.H_QSEQID]].strip()
        sacc = row[header_index[self.H_SSEQID]].strip()
        evalue_str = row[header_index[self.H_EVALUE]].strip()
        pident_str = row[header_index[self.H_PIDENT]].strip()
        qstart_str = row[header_index[self.H_QSTART]].strip()
        qend_str = row[header_index[self.H_QEND]].strip()
        sstart_str = row[header_index[self.H_SSTART]].strip()
        send_str = row[header_index[self.H_SEND]].strip()

        matched_seq = None
        if self.H_SSEQ in header_index and header_index[self.H_SSEQ] < len(row):
            cell = row[header_index[self.H_SSEQ]].strip()
            matched_seq = cell if cell != "" else None

        qstart = int(qstart_str)
        qend = int(qend_str)
        if qstart > qend:
            raise ValueError(f"qstart ({qstart}) must be <= qend ({qend}) in results TSV row: {row}")

        sstart = int(sstart_str)
        send = int(send_str)

        match = Match(
            query_accession=qacc,
            target_accession=sacc,
            query_start=qstart,
            query_end=qend,
            target_start=sstart,
            target_end=send,
            e_value=float(evalue_str.replace(",", "")) if evalue_str else None,
            identity=float(pident_str.replace(",", "")) if pident_str else None,
            matched_sequence=matched_seq,
        )
        return match


def extract_fragments(target_accession, dna_start, dna_end, seq):
    fragments = []
    current_start = 0

    seq += "*"
    for i in range(0, len(seq)):
        if seq[i] == "*":
            if i > current_start:
                start, end = to_dna_coordinate(dna_start, dna_end, current_start+1, i)
                fragments.append((target_accession, start, end, seq[current_start:i]))
            current_start = i+1

    return fragments


def get_aa_sequences(target_accession, target_sequence, target_left = None, target_right = None, strand = None, win=50000, win_overlap=10000, min_aa_length=8):

    # target_left and target_right are 1b, if provided
    if target_left is None:
        target_left = 1
    if target_right is None:
        target_right = len(target_sequence)

    coords = {}

    for win_i in range(target_left-1, target_right, win):
        # print(target_accession, target_left, target_right, "win", win_i+1, "+", win_overlap)

        subs = target_sequence[win_i:min(win_i+win+win_overlap, target_right)]
        translations = []
        if strand is None or strand == 1:
            translations.extend(compute_three_frame_translations(subs, 1, len(subs)))
        if strand is None or strand == -1:
            translations.extend(compute_three_frame_translations(subs, len(subs), 1))

        for translation in translations:
            frame_dna_start = translation[0]+win_i
            frame_dna_end = translation[1]+win_i

            for f in extract_fragments(target_accession, frame_dna_start, frame_dna_end, translation[2]):
                acc, start, end, seq = f
                k = (start, end)
                if len(seq) >= min_aa_length:
                    # print(k, seq)
                    coords[k] = f

    return list(coords.values())


def hmm_search_genome(hmm_file, genomic_fasta_dict, min_aa_length = 8,
                      target_accession = None, target_left = None, target_right = None,
                      strand = None, cpus = None):

    fragments = []
    for acc, genome_sequence in genomic_fasta_dict.items():
        if target_accession and acc != target_accession:
            continue
        fragments.extend(get_aa_sequences(acc, genome_sequence, target_left=target_left, target_right=target_right, strand=strand, min_aa_length=min_aa_length))

    fragments = [x for x in fragments if len(x[3]) >= min_aa_length]
    fragments = [x for x in fragments if (x[3].count('X') / len(x[3])) < 0.1]

    translated_fasta = {}
    name_to_coordinates = {}
    for target_accession, target_start, target_end, seq in fragments:
        target_name = f"{target_accession}_{target_start}_{target_end}"
        translated_fasta[target_name] = seq
        # print(target_name, seq)
        name_to_coordinates[target_name] = (target_accession, target_start, target_end)

    hmm_rows = hmmsearch_sequence_dict(hmm_file, translated_fasta, cpu=cpus)
    hmm_rows = [row for row in hmm_rows if row["dom_evalue"] <= DOM_EVALUE_LIMIT]

    filtered_rows = []
    for row in hmm_rows:
        # hmmsearch: query = hmm, target = protein
        hmm_name = row["query_accession"]
        if not hmm_name or hmm_name.strip() == "-":
            hmm_name = row["query_name"]
        target_name = row["target_name"]

        full_aa_seq = translated_fasta[target_name]
        aa_seq = extract_subsequence(full_aa_seq, row["ali_from"], row["ali_to"])
        if len(aa_seq) < min_aa_length:
            # print("throw away", len(aa_seq), "at target", row["ali_from"])
            continue

        target_accession, target_start, target_end = name_to_coordinates[target_name]
        dna_ali_from, dna_ali_to = to_dna_coordinate(target_start, target_end, row["ali_from"], row["ali_to"])

        row["query_accession"] = hmm_name
        row["target_accession"] = target_accession
        row["ali_from"] = dna_ali_from
        row["ali_to"] = dna_ali_to
        row["matched_sequence"] = aa_seq
        filtered_rows.append(row)

    return filtered_rows
