import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
from Bio.Seq import Seq

from .sequence import extract_subsequence, extract_subsequence_strand_sensitive, compute_three_frame_translations
from .match import Match, ProteinHit, order_matches_for_junctions, order_matches, NonlinearMatchException
from .detect import DOM_EVALUE_LIMIT, hmm_search_genome
from .hmm import hmmsearch, HMMCollection

from tangle.detected import DetectedTable
from tangle import unique_batch

VERBOSE_EXPAND = 0


@dataclass
class Candidate:
    assigned_overlap_to_left: Optional[int]  # for overlaps; None for gaps
    window_seq: str
    stitched: str
    left_trimmed: int
    right_kept: str


def generate_transition_candidates(
    left_aa: str,
    right_aa: str,
    overlap_len: int,
    gap_len: int,
    overlap_flanking_len: int = 40,
) -> List[Candidate]:

    candidates: List[Candidate] = []

    left_len = len(left_aa)
    right_len = len(right_aa)
    # have to add this here to deal with deletions in sequences - i.e. missing bps from query/profile
    overlap_len = min(left_len, right_len, overlap_len)

    assert left_len >= overlap_len
    assert right_len >= overlap_len
    assert overlap_len == 0 or gap_len == 0

    left_window_start = max(left_len-overlap_len-overlap_flanking_len, 0)
    right_window_end_plus1 = min(overlap_len+overlap_flanking_len, right_len)

    for k in range(0, overlap_len + 1):
        # Assign k of the overlap to the left, and (overlap_len - k) to the right.

        left_trim = (overlap_len - k)
        left_prefix = left_aa[: left_len - left_trim]
        left_window = left_aa[left_window_start: left_len - overlap_len + k]
        right_suffix = right_aa[k: right_len]
        right_window = right_aa[k: right_window_end_plus1]

        if overlap_len == 0 and gap_len:
            gap = "X" * gap_len
            assigned_overlap_to_left = None
        else:
            gap = ""
            assigned_overlap_to_left = k

        stitched = left_prefix + gap + right_suffix
        window_seq = left_window + gap + right_window

        candidates.append(
            Candidate(assigned_overlap_to_left=assigned_overlap_to_left,
                      window_seq=window_seq,
                      stitched=stitched, left_trimmed=left_trim, right_kept=gap+right_suffix))

    return candidates


def hmmsearch_find_best_candidate(hmm_file_name, sequences):
    matches = hmmsearch(hmm_file_name, sequences, cutoff=False, gap_removal=False)
    matches = [row for row in matches if row["seq_evalue"] <= DOM_EVALUE_LIMIT]

    best_idx = None
    best_score = float("-inf")
    best_evalue = None

    for match in matches:
        name = match["target_name"]
        score = match["seq_score"]
        evalue = match["seq_evalue"]

        if name.startswith("cand_"):
            idx = int(name.split("_")[1])
            if score > best_score:
                best_score = score
                best_evalue = evalue
                best_idx = idx

    return best_idx, best_score, best_evalue


def score_and_select_best_transition(
    candidates: List[Candidate],
    hmm_file_name: str,
) -> Candidate:

    if len(candidates) == 1:
        return candidates[0]
    sequences = [c.window_seq for c in candidates]
    best_idx, _1, _2 = hmmsearch_find_best_candidate(hmm_file_name, sequences)
    if best_idx is None:
        print("Cannot determine best candidate using hmmsearch, default to first transition")
        # print("candidates were", sequences)
        return candidates[0]
    return candidates[best_idx]


def aa_by_match(matches: List[Match]) -> Dict[int, str]:
    mapping: Dict[int, str] = {}
    for m in matches:
        aa_full = m.target_sequence_translated()
        mapping[id(m)] = aa_full
    return mapping


def stitch_cleaned_sequence(
    ordered_pairs: List[Tuple[Match, Match, int, int]],
    best_candidates_by_pair_index: Dict[int, Candidate],
    aa_by_match: Dict[int, str],
) -> str:

    result = ""
    for idx, (left, _right, _overlap, _gap) in enumerate(ordered_pairs):
        cand = best_candidates_by_pair_index[idx]
        if idx == 0:
            # print("stitch, add", cand.stitched)
            result += cand.stitched
        else:
            # print("stitch, trim", cand.left_trimmed)
            result = result[: len(result)-cand.left_trimmed]
            # print("stitch, add", cand.right_kept)
            result += cand.right_kept
    return result


def _clone(m: Match) -> Match:
    return Match(
        query_accession=m.query_accession,
        target_accession=m.target_accession,
        query_start=m.query_start,
        query_end=m.query_end,
        target_start=m.target_start,
        target_end=m.target_end,
        e_value=m.e_value,
        identity=m.identity,
        matched_sequence=m.matched_sequence,
        target_sequence=m.target_sequence,
    )


def _trim_dna_front(dna: Optional[str], aa_count: int) -> Optional[str]:
    if dna is None or aa_count <= 0:
        return dna
    bases = 3 * aa_count
    if bases >= len(dna):
        return ""
    return dna[bases:]


def _trim_dna_back(dna: Optional[str], aa_count: int) -> Optional[str]:
    if dna is None or aa_count <= 0:
        return dna
    bases = 3 * aa_count
    if bases >= len(dna):
        return ""
    return dna[: len(dna) - bases]


def adjust_target_coordinates(left: Match, right: Match, cand: Candidate) -> Tuple[Match, Match]:
    nl = _clone(left)
    nr = _clone(right)

    # Gap: leave as-is
    if cand.assigned_overlap_to_left is None:
        return nl, nr

    nl.query_end -= cand.left_trimmed
    nl.target_end -= 3 * cand.left_trimmed
    nl.target_sequence = _trim_dna_back(nl.target_sequence, cand.left_trimmed)

    nr.query_start += cand.assigned_overlap_to_left
    nr.target_start += 3 * cand.assigned_overlap_to_left
    nr.target_sequence = _trim_dna_front(nr.target_sequence, cand.assigned_overlap_to_left)

    # print("adjusting target coordinates, left", -cand.left_trimmed, "right", cand.assigned_overlap_to_left)
    # print("now left", nl.query_start, nl.query_end, "right", nr.query_start, nr.query_end)

    return nl, nr


def hmm_clean_protein(
    protein_hit: ProteinHit,
    hmm_file_name: str,
    overlap_flanking_len: int = 40,
    min_query_match_len: int = 4
) -> ProteinHit:

    if len(protein_hit.matches) < 2:
        new_protein_hit = ProteinHit(
            matches=protein_hit.matches,
            query_start=protein_hit.query_start,
            query_end=protein_hit.query_end,
            target_start=protein_hit.target_start,
            target_end=protein_hit.target_end,
            hmm_file=hmm_file_name
        )
        return new_protein_hit

    """
    print("")
    print("cleaning", protein_hit.protein_hit_id, protein_hit.query_accession)
    """

    old_matches = protein_hit.matches

    # remove very small AA matches
    old_matches = [m for m in old_matches if m.query_end-m.query_start+1 >= min_query_match_len]

    if len(old_matches) < 1:
        return None
    elif len(old_matches) == 1:
        cleaned_pm = ProteinHit(
            matches=old_matches,
            query_start=min(m.query_start for m in old_matches),
            query_end=max(m.query_end for m in old_matches),
            target_start=protein_hit.target_start,
            target_end=protein_hit.target_end,
            hmm_file=hmm_file_name
        )
        return cleaned_pm

    # Compute AA per match and junction candidates
    aa_map = aa_by_match(old_matches)
    pairs = order_matches_for_junctions(old_matches)

    selected: Dict[int, Candidate] = {}
    for idx, (left, right, overlap_len, gap_len) in enumerate(pairs):
        """
        print(idx, "left", left.query_start, left.query_end, left.target_start, left.target_end, aa_map[id(left)])
        print(idx, "right", right.query_start, right.query_end, right.target_start, right.target_end, aa_map[id(right)])
        print(idx, "overlap/gap", overlap_len, gap_len)
        """
        cands = generate_transition_candidates(
            aa_map[id(left)], aa_map[id(right)], overlap_len, gap_len, overlap_flanking_len
        )
        # print("choosing candidate for", left.query_start, left.query_end, "and", right.query_start, right.query_end)
        best = cands[0] if len(cands) <= 1 else score_and_select_best_transition(cands, hmm_file_name)
        selected[idx] = best
        # print("  chose", best)

    # Stitch the final AA from original matches and chosen splits
    cleaned_aa = stitch_cleaned_sequence(pairs, selected, aa_map)

    # Create new Match objects
    new_matches: List[Match] = []
    current_left = pairs[0][0]
    for idx, (_left, right, _1, _2) in enumerate(pairs):
        selected_candidate = selected[idx]
        new_left, new_right = adjust_target_coordinates(current_left, right, selected_candidate)
        new_matches.append(new_left)
        current_left = new_right
    new_matches.append(current_left)

    if not ProteinHit.can_collate_from_matches(new_matches):
        print("Cleaned matches cannot be collated, revert")
        print(protein_hit.collated_protein_sequence)
        return protein_hit

    cleaned_pm = ProteinHit(
        matches=new_matches,
        query_start=min(m.query_start for m in new_matches),
        query_end=max(m.query_end for m in new_matches),
        target_start=protein_hit.target_start,
        target_end=protein_hit.target_end,
        hmm_file=hmm_file_name
    )
    # Validate cleaned sequence matches the newly collated sequence from adjusted matches
    if cleaned_pm.collated_protein_sequence != cleaned_aa:
        print(f"Cleaned ProteinHit collated sequence mismatch: {cleaned_pm.collated_protein_sequence} != {cleaned_aa}")
    return cleaned_pm


def hmm_clean(protein_hits: List[ProteinHit], hmm_collection: HMMCollection, overlap_flanking_len: int = 40) -> List[ProteinHit]:

    cleaned: Dict[ProteinHit] = {}

    skipped = []
    for pm in protein_hits:
        hmm_profile = hmm_collection.get(pm.query_accession)
        if hmm_profile is None:
            if pm.query_accession not in skipped:
                skipped.append(pm.query_accession)
                print("skipping", pm.query_accession, "cannot find HMM profile")
        else:
            new_pm = hmm_clean_protein(pm, hmm_profile, overlap_flanking_len)
            if new_pm:
                cleaned[new_pm.protein_hit_id] = new_pm

    return list(cleaned.values())


def find_more_matches_at_locus(query_accession, hmm_file, old_start, old_end, target_full_sequence, target_accession, target_left, target_right, strand, cpus=None):

    hmm_rows = hmm_search_genome(
        hmm_file, {target_accession: target_full_sequence},
        target_accession = target_accession,
        target_left = target_left,
        target_right = target_right,
        strand = strand,
        cpus = cpus,
        filter_by_evalue_cond = True  # we already assume there's a protein here...
    )

    if not hmm_rows:
        if VERBOSE_EXPAND:
            print("no hmmsearch results")
        return None

    new_matches = []
    for hmm_row in hmm_rows:
        target_sequence = extract_subsequence_strand_sensitive(target_full_sequence, hmm_row["ali_from"], hmm_row["ali_to"])
        match = Match(
            query_accession=query_accession,
            target_accession=target_accession,
            query_start=hmm_row["hmm_from"],
            query_end=hmm_row["hmm_to"],
            target_start=hmm_row["ali_from"],
            target_end=hmm_row["ali_to"],
            e_value=hmm_row["dom_evalue"],
            identity=None,
            target_sequence=target_sequence,
            matched_sequence=hmm_row["matched_sequence"]
        )
        assert match.matched_sequence == match.target_sequence_translated()
        new_matches.append(match)

    index_of_old_start = None
    index_of_old_end = None
    index_of_new_start = None
    index_of_new_end = None

    if strand > 0:
        new_matches = sorted(new_matches, key=lambda m: m.target_start)
        for i, match in enumerate(new_matches):
            if match.target_end > old_start and index_of_old_start is None:
                index_of_old_start = i
            if match.target_end >= old_end and index_of_old_end is None:
                index_of_old_end = i
        if index_of_old_end is None and index_of_old_start is not None:
            index_of_old_end = len(new_matches)-1

        if index_of_old_start is None or index_of_old_end is None or index_of_old_end < index_of_old_start:
            if VERBOSE_EXPAND:
                print("cannot find old start and end indices on fwd strand")
                for i, match in enumerate(new_matches):
                    mark = "   "
                    if i in (index_of_old_start, index_of_old_end):
                        mark = " ->"
                    if i in (index_of_new_start, index_of_new_end):
                        mark = " =>"
                    print(f" {mark} {match.target_start}, {match.target_end}, {match.query_start}, {match.query_end}")
            return None

    else:
        new_matches = sorted(new_matches, key=lambda m: -m.target_start)
        for i, match in enumerate(new_matches):
            if match.target_end < old_start and index_of_old_start is None:
                index_of_old_start = i
            elif match.target_end <= old_end and index_of_old_end is None:
                index_of_old_end = i
        if index_of_old_end is None and index_of_old_start is not None:
            index_of_old_end = len(new_matches)-1

        if index_of_old_start is None or index_of_old_end is None or index_of_old_end < index_of_old_start:
            if VERBOSE_EXPAND:
                print("cannot find old start and end indices on rev strand")
                for i, match in enumerate(new_matches):
                    mark = "   "
                    if i in (index_of_old_start, index_of_old_end):
                        mark = " ->"
                    if i in (index_of_new_start, index_of_new_end):
                        mark = " =>"
                    print(f" {mark} {match.target_start}, {match.target_end}, {match.query_start}, {match.query_end}")
            return None
  
    # sanity check - if we cannot collate what we think are the old matches,
    # that means what we think are the old matches are likely repeats of the
    # same match so we can just stop

    if not ProteinHit.can_collate_from_matches(new_matches[index_of_old_start:index_of_old_end+1]):
        if VERBOSE_EXPAND:
            print("sticking with the old matches")
        return None
 
    # move starting index as far back as we can
    index_of_new_start = index_of_old_start
    i = index_of_old_start-1
    while i >= 0:
        if ProteinHit.can_collate_from_matches(new_matches[i:index_of_old_end+1]):
            index_of_new_start = i
            i -= 1
        else:
            break
    
    # move ending index as far forward as we can
    index_of_new_end = index_of_old_end
    i = index_of_old_end+1
    while i <= len(new_matches)-1:
        if ProteinHit.can_collate_from_matches(new_matches[index_of_old_start:i+1]):
            index_of_new_end = i
            i += 1
        else:
            break

    if VERBOSE_EXPAND:
        print("found:")
        for i, match in enumerate(new_matches):
            mark = "  "
            if i in (index_of_old_start, index_of_old_end):
                mark = "->"
            if i in (index_of_new_start, index_of_new_end):
                mark = "=>"
            print(f"  {mark} {match.target_start}, {match.target_end}, {match.query_start}, {match.query_end}")

    return new_matches[index_of_new_start:index_of_new_end+1]


def hmm_expand_protein(protein_hit, genomic_sequence_dict, hmm_file, cpus = None):

    target_full_sequence = genomic_sequence_dict[protein_hit.target_accession]
    query_accession = protein_hit.matches[0].query_accession
    target_accession = protein_hit.matches[0].target_accession
    strand = -1 if protein_hit.matches[0].on_reverse_strand else 1
    start = protein_hit.target_start
    end = protein_hit.target_end

    max_search_distance = 30000
    target_left = max(min(start, end) - max_search_distance, 1)
    target_right = min(max(start, end) + max_search_distance, len(target_full_sequence))

    if VERBOSE_EXPAND:
        print(f"{query_accession} on {target_accession}, {target_left}-{target_right} (based on {start}-{end}), strand {strand}, contig {len(target_full_sequence)}")
        for i, match in enumerate(protein_hit.matches):
            print(f" old {match.target_start}, {match.target_end}, {match.query_start}, {match.query_end}")

    new_matches = find_more_matches_at_locus(
        query_accession, hmm_file, start, end,
        target_full_sequence, target_accession, target_left, target_right, strand,
        cpus = cpus
    )

    if new_matches is None:
        return protein_hit

    new_pm = ProteinHit(
        matches=new_matches,
        query_start=min(m.query_start for m in new_matches),
        query_end=max(m.query_end for m in new_matches),
        target_start=min(m.target_start for m in new_matches) if protein_hit.target_start < protein_hit.target_end else max(m.target_start for m in new_matches),
        target_end=max(m.target_end for m in new_matches) if protein_hit.target_start < protein_hit.target_end else min(m.target_end for m in new_matches),
        hmm_file=hmm_file
    )

    return new_pm


def hmm_expand(protein_hits, genomic_sequence_dict, hmm_collection, cpus = None):
    new_protein_hits = {}

    skipped = []
    for pm in protein_hits:
        """
        print()
        print(pm.protein_hit_id)
        for nm in sorted(pm.matches, key=lambda m: m.target_start):
            print("    ", nm.target_start, nm.target_end, nm.query_start, nm.query_end)
        """

        hmm_profile = hmm_collection.get(pm.query_accession)
        if hmm_profile is None:
            if pm.query_accession not in skipped:
                skipped.append(pm.query_accession)
                print("skipping", pm.query_accession, "cannot find HMM profile")
        else:
            pm = hmm_expand_protein(pm, genomic_sequence_dict, hmm_profile, cpus = cpus)
            new_protein_hits[pm.protein_hit_id] = pm

    return list(new_protein_hits.values())


def write_fasta_record(f, pm: ProteinHit) -> None:
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
    genome_accession: str,
    protein_hits: List[ProteinHit],
    proteins_aa_path: str,
    proteins_tsv_path: str,
    append = False
) -> None:

    protein_hits = [pm for pm in protein_hits if pm.can_produce_single_sequence()]

    with open(proteins_aa_path, "at" if append else "wt") as f_prot:
        for pm in protein_hits:
            write_fasta_record(f_prot, pm)

    write_tsv(proteins_tsv_path, protein_hits, genome_accession, append)
