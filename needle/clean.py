import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from .match import Match, ProteinHit, order_matches_for_junctions, NonlinearMatchException
from .detect import DOM_EVALUE_LIMIT
from .hmm import hmmsearch, HMMCollection


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
    try:
        pairs = order_matches_for_junctions(old_matches)
    except NonlinearMatchException as e:
        print("Cannot order existing detected matches, revert")
        print(protein_hit.collated_protein_sequence)
        return protein_hit

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
