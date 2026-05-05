import os
import csv
import errno
import hashlib
import itertools
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
from Bio.Seq import Seq
import networkx as nx

MAX_AA_OVERLAP_GROUPING = 15


@dataclass(frozen=True)
class Match:  # does not support matches across circular boundary
    query_accession: str
    target_accession: str
    query_start: int
    query_end: int
    target_start: int  # 1-based, 5' to 3', so target_start > target_end on reverse strand
    target_end: int    # 1-based, 5' to 3', so target_start > target_end on reverse strand
    e_value: float
    matched_sequence: Optional[str] = None
    target_sequence: Optional[str] = None  # 5' to 3'
    score: Optional[float] = None

    @property
    def on_reverse_strand(self) -> bool:
        return self.target_start > self.target_end

    @property
    def sortable_target_start(self):
        return -self.target_start if self.on_reverse_strand else self.target_start

    @property
    def sortable_target_end(self):
        return -self.target_end if self.on_reverse_strand else self.target_end

    def query_len(self):
        return abs(self.query_end - self.query_start) + 1

    def target_sequence_translated(self) -> str:
        if not self.target_sequence or self.query_start > self.query_end:
            return ""
        dna = self.target_sequence.upper()
        usable_len = (len(dna) // 3) * 3
        if usable_len == 0:
            return ""
        dna = dna[:usable_len]
        return str(Seq(dna).translate(table="Standard", to_stop=False))


class NonlinearMatchException(Exception):
    pass


def remove_duplicate_matches(matches: List[Match]) -> List[Match]:
    keyf = lambda m: (m.query_accession, m.query_start, m.query_end, m.target_accession, m.target_start, m.target_end)
    matches = sorted(matches, key=keyf)
    retained = []
    for k,group in itertools.groupby(matches, keyf):
        retained.append(list(group)[0])
    return retained


class Reason(object):

    def __init__(self):
        self.reason = None

    def was(self, msg):
        self.reason = msg


def match_is_pred_of(left, right, reason):
    """
    returns True if left is a pred of right. requires left and right being
    ordered by sortable_target_start, which is equivalent to saying they really
    should be ordered by query coordinates.
    """

    assert left.sortable_target_start <= right.sortable_target_start

    # matches should be on same strand
    if left.on_reverse_strand != right.on_reverse_strand:
        reason.was("Matches are on different strands")
        return False

    if left.query_start > right.target_start:
        msg = "Consecutive query matches are reversed on DNA, new copy? (%s, %s, rev strand %s, dna %s-%s (hmm %s-%s %s) and %s-%s (hmm %s-%s %s))" % \
                (left.query_accession, left.target_accession, left.on_reverse_strand,
                 left.target_start, left.target_end, left.query_start, left.query_end, left.e_value,
                 right.target_start, right.target_end, right.query_start, right.query_end, right.e_value)
        reason.was(msg)
        return False

    # query cannot contain each other
    if right.query_end <= left.query_end:
        msg = "Matching queries contain each other (%s %s-%s and %s-%s)" % \
                (left.target_accession, left.query_start, left.query_end, right.query_start, right.query_end)
        reason.was(msg)
        return False

    # left cannot contain right on target
    if (left.on_reverse_strand is False and right.target_end < left.target_end) or \
       (left.on_reverse_strand is True and right.target_end > left.target_end):
        msg = "DNA matches contain each other"
        reason.was(msg)
        return False

    # query overlap not bigger than either left or right sequence

    query_overlap_len = max(0, left.query_end - right.query_start + 1)
    if query_overlap_len > len(left.target_sequence_translated()) or \
       query_overlap_len > len(right.target_sequence_translated()):
        msg = "Overlap larger than matched sequence"
        reason.was(msg)
        return False

    return True


def build_graph(matches, reasons = None):

    graph = nx.DiGraph()
    matches = sorted(matches, key=lambda m: (m.sortable_target_start, m.sortable_target_end))

    prev_matches = []

    for match in matches:
        found_pred = False

        for possible_pred in prev_matches[::-1]:
            # already linked
            if graph.has_node(possible_pred) and \
               graph.has_node(match) and \
               nx.has_path(graph, possible_pred, match):
                continue

            reason = Reason()
            if match_is_pred_of(possible_pred, match, reason):
                ancestors = nx.ancestors(graph, possible_pred)
                if len(ancestors) == 0 or \
                   match.query_start > max([x.query_end for x in ancestors]):
                    graph.add_edge(possible_pred, match)
                    found_pred = True
                elif reasons is not None:
                   reason.was("Junctions overlap")
                   reasons.append(reason)
            elif reasons is not None:
                reasons.append(reason)

        if not found_pred:
            graph.add_node(match)

        prev_matches.append(match)

    return graph


def graph_paths(graph):

    roots = [n for n, d in graph.in_degree() if d == 0]
    leaves = [n for n, d in graph.out_degree() if d == 0]

    all_paths = []
    for root in roots:
        paths = nx.all_simple_paths(graph, source=root, target=leaves)
        all_paths.extend(paths)

    return all_paths


def order_matches(matches: List[Match]) -> List[Match]:
    if not matches:
        return []

    reasons = []
    graph = build_graph(matches, reasons)
    paths = graph_paths(graph)
    if len(paths) != 1:
        msg = "; ".join([r.reason for r in reasons])
        raise NonlinearMatchException(msg)

    return sorted(matches, key=lambda m: (m.query_start, m.query_end))


def order_matches_for_junctions(matches: List[Match]) -> List[Tuple[Match, Match, int, int]]:
    if not matches:
        return []

    ordered = order_matches(matches)
    pairs: List[Tuple[Match, Match, int, int]] = []
    junctions: List[Tuple[int, int]] = []

    i = 0
    while i < len(ordered)-1:
        left = ordered[i]
        right = ordered[i + 1]
        assert right.query_start >= left.query_start  # sort condition

        query_overlap_len = max(0, left.query_end - right.query_start + 1)
        gap_len = max(0, right.query_start - left.query_end - 1)
        pairs.append((left, right, query_overlap_len, gap_len))
        if gap_len:
            junctions.append((left.query_end, right.query_start))
        else:
            junctions.append((right.query_start, left.query_end))
        i += 1

    return pairs


@dataclass
class ProteinHit:
    matches: List[Match]
    query_start: int
    query_end: int
    target_start: int  # 1-based, 5' to 3' of gene, so target_start > target_end on reverse strand
    target_end: int    # 1-based, 5' to 3' of gene, so target_start > target_end on reverse strand
    hmm_cleaned_protein_sequence: Optional[str] = None
    hmm_file: Optional[str] = None
    _protein_hit_id: Optional[str] = None
    _product_name: Optional[str] = None

    @property
    def on_reverse_strand(self) -> bool:
        return self.target_start > self.target_end

    @staticmethod
    def can_collate_from_matches(matches, verbose=False) -> bool:
        try:
            pairs = order_matches_for_junctions(matches)
        except NonlinearMatchException as e:
            if verbose:
                print(str(e))
            return False
        return True

    def can_collate(self) -> bool:
        return self.can_collate_from_matches(self.matches)

    @staticmethod
    def can_produce_single_sequence_from_matches(matches) -> bool:
        try:
            pairs = order_matches_for_junctions(matches)
        except NonlinearMatchException:
            return False
        for left, right, overlap, gaps in pairs:
            if overlap > 0:
                return False
        return True

    def can_produce_single_sequence(self) -> bool:
        return self.can_produce_single_sequence_from_matches(self.matches)

    @property
    def collated_protein_sequence(self) -> str:
        collated = ""
        if len(self.matches) < 2:
            return self.matches[0].target_sequence_translated()
        pairs = order_matches_for_junctions(self.matches)
        cur_left_aa = pairs[0][0].target_sequence_translated()
        for left, right, overlap, gaps in pairs:
            right_aa = right.target_sequence_translated()
            if gaps:
                new_s = cur_left_aa
                new_s += "X" * gaps
                collated += new_s
                cur_left_aa = right_aa
            else:
                new_s = cur_left_aa[0:len(cur_left_aa)-overlap]
                if overlap > 0:
                    new_s += "("+cur_left_aa[len(cur_left_aa)-overlap:]+"/"+right_aa[0:overlap]+")"
                collated += new_s
                cur_left_aa = right_aa[overlap:]
        collated += cur_left_aa
        return collated

    @property
    def protein_hit_id(self) -> str:
        if self._protein_hit_id is not None:
            return self._protein_hit_id
        assert self.matches
        ordered = sorted(self.matches, key=lambda m: (m.query_start, m.query_end))
        first = ordered[0]
        base_q = first.query_accession
        base_t = first.target_accession
        hasher = hashlib.sha1()
        for m in ordered:
            parts = [
                m.query_accession,
                m.target_accession,
                str(m.query_start),
                str(m.query_end),
                str(m.target_start),
                str(m.target_end),
            ]
            hasher.update("|".join(parts).encode("utf-8"))
        digest8 = hasher.hexdigest()[:8]
        # target accession which is contig accession may have multiplexed "|"
        # in it to be stripped off, so we put that at the beginning
        self._protein_hit_id = f"{base_t}_{base_q}_{digest8}"
        return self._protein_hit_id

    @property
    def query_accession(self) -> str:
        assert self.matches
        first = sorted(self.matches, key=lambda m: (m.query_start, m.query_end))[0]
        return first.query_accession

    @property
    def target_accession(self) -> str:
        assert self.matches
        first = sorted(self.matches, key=lambda m: (m.query_start, m.query_end))[0]
        return first.target_accession


def group_matches(all_matches, max_intron_length: int = 10_000, max_overlap_len: int = MAX_AA_OVERLAP_GROUPING) -> List[List[Match]]:
    """
    cluster Match objects
      - groups by (query_accession, target_accession)
      - within a (query, target) group, further splits into clusters if adjacent
        matches on the target are separated by more than max_intron_length.
    """

    if not all_matches:
        return []

    def target_interval(m: Match) -> (int, int):
        return (min(m.target_start, m.target_end), max(m.target_start, m.target_end))

    grouped: Dict[tuple, List[Match]] = {}
    for m in all_matches:
        key = (m.query_accession, m.target_accession, m.on_reverse_strand)
        grouped.setdefault(key, []).append(m)

    final_clusters = []

    for (query_id, target_id, on_reverse_strand), group in grouped.items():
        # Sort by target interval start to create distance-based clusters
        group_sorted_by_target = sorted(group, key=lambda m: target_interval(m)[0])

        clusters: List[List[Match]] = []
        current_cluster: List[Match] = []
        current_right: Optional[int] = None
        current_query = None

        # print("grouping", query_id, "on", target_id, "rev", on_reverse_strand)
        for m in group_sorted_by_target:
            left_t, right_t = target_interval(m)
            fragment_len = m.query_end - m.query_start + 1

            # print("  left", left_t, "right", right_t, "current query", current_query, "match", m.query_start, m.query_end)

            # New cluster
            if not current_cluster:
                current_cluster = [m]
                current_right = right_t
                current_query = (m.query_start, m.query_end)
                # print("    start new cluster")

            else:
                distance = left_t - (current_right or left_t)
                # Too far by target nuc distance, start new cluster
                if distance > max_intron_length:
                    clusters.append(current_cluster)
                    current_cluster = [m]
                    current_right = right_t
                    current_query = (m.query_start, m.query_end)
                    # print("    too far on target, start new cluster")

                elif distance < 0:  # overlap on the genome, we can't distinguish the overlapping copies so just group them together for now
                    current_cluster.append(m)
                    current_right = max(current_right or right_t, right_t)
                    current_query = (m.query_start, m.query_end)
                    # print("    overlap on genome, add to current cluster")

                elif current_right > right_t:  # this one is contained in the last match, start new cluster
                    clusters.append(current_cluster)
                    current_cluster = [m]
                    current_right = right_t
                    current_query = (m.query_start, m.query_end)
                    # print("    contained match? start new cluster")

                # Query rewound, start new cluster
                # Criteria here is:
                #    There is overlap on query and not overlap on target, and
                #      Too long overlap or
                #      Overlap is >50% of fragment
                #      Completely rewond (e.g. repeat)
                elif (on_reverse_strand is False and m.query_start < current_query[1] and \
                      ((current_query[1] - m.query_start + 1) > max_overlap_len or \
                       (current_query[1] - m.query_start + 1) / fragment_len > 0.5 or \
                       m.query_start <= current_query[0])) \
                    or \
                     (on_reverse_strand is True and m.query_end > current_query[0] and \
                      ((m.query_end - current_query[0] + 1) > max_overlap_len or \
                       (m.query_end - current_query[0] + 1) / fragment_len > 0.5 or \
                       m.query_start >= current_query[0])):
                    clusters.append(current_cluster)
                    current_cluster = [m]
                    current_query = (m.query_start, m.query_end)
                    # print("    query rewound, start new cluster")

                # Add to cluster
                else:
                    current_cluster.append(m)
                    current_right = max(current_right or right_t, right_t)
                    current_query = (m.query_start, m.query_end)
                    # print("    add to current cluster")

        if current_cluster:
            clusters.append(current_cluster)

        for cluster in clusters:
            cluster = remove_duplicate_matches(cluster)
            final_clusters.append(cluster)

    return final_clusters
