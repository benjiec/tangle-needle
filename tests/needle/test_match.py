import os
import tempfile
import unittest
import random

from needle.match import group_matches, ProteinHit, Match, NonlinearMatchException
from needle.match import Reason, match_is_pred_of, build_graph, graph_paths, order_matches_by_query, ordered_pairs


class TesterBase(unittest.TestCase):

    @staticmethod
    def makeM(query_start, query_end, target_start, target_end, query_accession=None, target_accession=None, target_sequence=None):
        if target_sequence is None:
            target_sequence="A"*(abs(target_end-target_start)+1)

        return Match(
            query_accession=query_accession, target_accession=target_accession, e_value=0,
            query_start=query_start, query_end=query_end, target_start=target_start, target_end=target_end,
            target_sequence=target_sequence
        )


class TestIsPred(TesterBase):

    def test_match_is_pred_of(self):
        m1 = self.makeM(1, 10, 1, 30)
        m2 = self.makeM(8, 15, 24, 45)
        m3 = self.makeM(18, 20, 60, 66)

        r = Reason()
        self.assertEqual(match_is_pred_of(m1, m2, r), True)
        self.assertEqual(match_is_pred_of(m2, m3, r), True)
        self.assertEqual(match_is_pred_of(m1, m3, r), True)

    def test_not_pred_if_reversed_on_query(self):
        m1 = self.makeM(8, 19, 1, 30)
        m2 = self.makeM(1, 10, 24, 45)

        r = Reason()
        self.assertEqual(match_is_pred_of(m1, m2, r), False)
        self.assertEqual(r.reason.startswith("Matched to query in unexpected direction"), True)

        m1 = self.makeM(8, 19, 60, 31)
        m2 = self.makeM(1, 10, 35, 1)

        r = Reason()
        self.assertEqual(match_is_pred_of(m1, m2, r), False)
        self.assertEqual(r.reason.startswith("Matched to query in unexpected direction"), True)

    def test_not_pred_if_query_regions_contained_in_one(self):
        # aaaaaaaaa
        #      bbb

        m1 = self.makeM(1, 10, 1, 30)
        m2 = self.makeM(6, 9, 24, 45)

        r = Reason()
        self.assertEqual(match_is_pred_of(m1, m2, r), False)
        self.assertEqual(r.reason.startswith("One matched query region contains the other"), True)

    def test_not_pred_if_target_regions_contained_in_one(self):
        # aaaaaaaaa
        #      bbb

        m1 = self.makeM(1, 10, 30, 60)
        m2 = self.makeM(6, 11, 33, 65)

        r = Reason()
        self.assertEqual(match_is_pred_of(m1, m2, r), True)

        m2 = self.makeM(6, 11, 33, 45)
        self.assertEqual(match_is_pred_of(m1, m2, r), False)
        self.assertEqual(r.reason.startswith("One matched target region contains the other"), True)

        # now reversed strand

        m1 = self.makeM(1, 10, 60, 30)
        m2 = self.makeM(6, 11, 45, 25)

        self.assertEqual(match_is_pred_of(m1, m2, r), True)

        m2 = self.makeM(6, 11, 45, 33)
        self.assertEqual(match_is_pred_of(m1, m2, r), False)
        self.assertEqual(r.reason.startswith("One matched target region contains the other"), True)

    def test_not_pred_if_query_coordinates_do_not_align_with_target_coordinates(self):

        m1 = self.makeM(1, 10, 1, 30)
        m2 = self.makeM(6, 16, 48, 30)

        r = Reason()
        self.assertEqual(match_is_pred_of(m2, m1, r), False)
        self.assertEqual(r.reason, "Matches are on different strands")

    def test_not_pred_if_query_overlap_is_larger_than_left_or_right_sequence(self):

        m1 = self.makeM(6,  16, 1,  21)
        m2 = self.makeM(10, 20, 25, 50)

        r = Reason()
        # overlap is 7 aa, so nuc sequence should be 21 bp or longer
        self.assertEqual(match_is_pred_of(m1, m2, r), True)

        # make nc sequence 20 bps, so overlap is longer, not allowed
        m1 = self.makeM(6,  16, 1,  20)
        self.assertEqual(match_is_pred_of(m1, m2, r), False)
        self.assertEqual(r.reason.startswith("Overlap between matched sequences longer than one of the matched sequence"), True)

    def test_not_pred_if_overlap_of_matched_query_regions_is_too_large(self):

        m1 = self.makeM( 5, 35, 1001, 1090)
        m2 = self.makeM(35, 65, 1101, 1190)
        m3 = self.makeM(65, 95, 1201, 1290)
        m4 = self.makeM(70, 99, 1301, 1390)

        r = Reason()
        self.assertEqual(match_is_pred_of(m1, m2, r, max_aa_overlap = 20), True)
        self.assertEqual(match_is_pred_of(m2, m3, r, max_aa_overlap = 20), True)
        self.assertEqual(match_is_pred_of(m3, m4, r, max_aa_overlap = 20), False)
        self.assertEqual(r.reason.startswith("Overlap between matches longer than threshold 20 aa"), True)

    def test_not_pred_if_matched_query_regions_overlap_completely(self):

        m1 = self.makeM( 5, 18, 1001, 1020)
        m2 = self.makeM(15, 28, 1101, 1120)
        m3 = self.makeM(25, 38, 1201, 1220)
        m4 = self.makeM(25, 48, 1301, 1520)

        r = Reason()
        self.assertEqual(match_is_pred_of(m1, m2, r), True)
        self.assertEqual(match_is_pred_of(m2, m3, r), True)
        self.assertEqual(match_is_pred_of(m3, m4, r), False)
        self.assertEqual(r.reason.startswith("One matched query region contains the other"), True)


class TestGraphPaths(TesterBase):

    def test_produces_single_path_if_matches_are_successors(self):
        m1 = self.makeM(1, 10, 1, 30)
        m2 = self.makeM(8, 15, 24, 45)
        m3 = self.makeM(18, 20, 60, 66)

        graph = build_graph([m1, m2, m3])
        paths = graph_paths(graph)
        self.assertEqual(len(paths), 1)
        self.assertEqual(paths[0], [m1, m2, m3])

    def test_produces_multiple_paths_if_matches_are_not_successors(self):
        # aaaaaaaaa
        #      bbbbbbbbbbb
        #            ccc
        #                ddd

        m1 = self.makeM(1, 10, 1, 30)
        m2 = self.makeM(6, 16, 16, 48)
        m3 = self.makeM(12, 14, 34, 42)
        m4 = self.makeM(16, 18, 46, 54)

        graph = build_graph([m1, m2, m3, m4])
        paths = graph_paths(graph)
        self.assertEqual(len(paths), 2)
        self.assertCountEqual(
            paths,
            [[m1, m2, m4], [m1, m3, m4]]
        )

    def test_produces_multiple_paths_if_junctions_overlap(self):
        # aaaaaaaa
        #      bbbbbb
        #        cccccccc

        m1 = self.makeM(1, 8,  1, 24)
        m2 = self.makeM(6, 12, 16, 36)
        m3 = self.makeM(8, 16, 22, 48)

        graph = build_graph([m1, m2, m3])
        paths = graph_paths(graph)
        self.assertEqual(len(paths), 2)
        self.assertCountEqual(
            paths,
            [[m1, m2], [m1, m3]]
        )

    def test_multiple_starting_nodes_for_paths(self):
        # aaaaaaaaa
        #   bbb ccc
        #           ddd

        m1 = self.makeM(1,  20, 1,   60)
        m2 = self.makeM(6,  10, 70,  80)
        m3 = self.makeM(12, 20, 90,  104)
        m4 = self.makeM(21, 24, 120, 130)

        graph = build_graph([m1, m2, m3, m4])
        paths = graph_paths(graph)
        self.assertEqual(len(paths), 2)
        self.assertCountEqual(
            paths,
            [[m1, m4], [m2, m3, m4]]
        )

    def test_does_not_connect_node_if_possible_pred_already_has_descendants_covering_similar_region(self):
        # aaaaaaaaa
        #           bbb
        #                      ccc

        m1 = self.makeM(1,  20, 1,   60)
        m2 = self.makeM(20, 30, 70,  100)
        m3 = self.makeM(20, 30, 100, 130)

        graph = build_graph([m1, m2, m3])
        paths = graph_paths(graph)
        self.assertEqual(len(paths), 2)
        self.assertCountEqual(
            paths,
            [[m1, m2], [m3]]
        )


class TestOrderingMatches(TesterBase):

    def test_order_matches_by_query_errors_if_multiple_paths(self):
        m1 = self.makeM(1,  10, 1,  30)
        m2 = self.makeM(8,  15, 24, 45)
        m3 = self.makeM(18, 20, 60, 66)

        r = order_matches_by_query([m1, m3, m2])  # input not in order
        self.assertEqual(r, [m1, m2, m3])

        m1 = self.makeM(1,  10, 1,  30)
        m2 = self.makeM(8,  9,  24, 45)
        m3 = self.makeM(18, 20, 60, 66)

        with self.assertRaisesRegex(NonlinearMatchException, "One matched query region contains the other"):
            r = order_matches_by_query([m1, m3, m2])  # input not in order

    def test_ordered_pairs_return_junction_information(self):
        m1 = self.makeM(1,  10, 1,  30)
        m2 = self.makeM(8,  15, 24, 45)
        m3 = self.makeM(18, 20, 60, 66)

        pairs = ordered_pairs([m1, m3, m2])  # input not in order
        self.assertEqual(len(pairs), 2)
        self.assertEqual(pairs[0], (m1, m2, 3, 0))
        self.assertEqual(pairs[1], (m2, m3, 0, 2))


class TestGroupingMatches(TesterBase):

    def test_group_matches_separate_matches_by_contig_and_distance(self):

        m1 = self.makeM( 5, 10, 1001, 1020, "Q", "S1")
        m2 = self.makeM(15, 20, 1101, 1120, "Q", "S1")
        m3 = self.makeM(25, 30, 1201, 1220, "Q", "S1")
        m4 = self.makeM(35, 40, 1301, 1320, "Q", "S2")
        m5 = self.makeM(45, 50, 1401, 1420, "Q", "S2")
        m6 = self.makeM(55, 60, 1501, 1520, "Q", "S2")
        m7 = self.makeM(65, 70, 11601, 11620, "Q", "S2")
        m8 = self.makeM(70, 80, 11701, 11720, "Q", "S2")
        m9 = self.makeM(75, 90, 11801, 11820, "Q", "S2")

        matches = [m1, m2, m3, m4, m5, m6, m7, m8, m9]
        random.shuffle(matches)
        clusters = group_matches(matches)

        self.assertCountEqual(clusters, [
            [m1, m2, m3],
            [m4, m5, m6],
            [m7, m8, m9]
        ])

    def test_group_matches_separate_matches_on_reverse_strand(self):

        m1 = self.makeM(25, 30, 1020, 1001, "Q", "S1")
        m2 = self.makeM(15, 20, 1120, 1101, "Q", "S1")
        m3 = self.makeM( 5, 10, 1220, 1201, "Q", "S1")
        m4 = self.makeM(35, 40, 11301, 11320, "Q", "S1")
        m5 = self.makeM(45, 50, 11401, 11420, "Q", "S1")
        m6 = self.makeM(55, 60, 11501, 11520, "Q", "S1")

        matches = [m1, m2, m3, m4, m5, m6]
        random.shuffle(matches)
        clusters = group_matches(matches)

        self.assertCountEqual(clusters, [
            [m3, m2, m1],
            [m4, m5, m6],
        ])

    def test_group_matches_separate_matches_on_different_strands(self):
        m1 = self.makeM( 5, 18, 1001, 1020, "Q", "S1")
        m2 = self.makeM(15, 28, 1101, 1120, "Q", "S1")
        m3 = self.makeM(25, 38, 1201, 1220, "Q", "S1")

        matches = [m1, m2, m3]
        random.shuffle(matches)
        clusters = group_matches(matches)
        self.assertCountEqual(clusters, [
            [m1, m2, m3],
        ])
        
        m3 = self.makeM(25, 38, 1220, 1201, "Q", "S1")
        matches = [m1, m2, m3]
        random.shuffle(matches)
        clusters = group_matches(matches)
        self.assertCountEqual(clusters, [
            [m1, m2],
            [m3]
        ])


class TestProteinHitClass(TesterBase):

    def test_protein_hit_id_deterministic_and_changes_when_inputs_change(self):
        # Two identical ProteinHit objects should have the same ID
        a1 = Match("QID","TID",1,3,1,9,0.0,100.0)
        a2 = Match("QID","TID",4,6,10,18,1e-12,95.0)
        pm1 = ProteinHit([a1, a2], 1, 6, 1, 18)
        pm2 = ProteinHit([a1, a2], 1, 6, 1, 18)
        self.assertEqual(pm1.protein_hit_id, pm2.protein_hit_id)
        # Change an input value (e.g., query_start) to produce a different ID
        b2 = Match("QID","TID",5,6,10,18,1e-12,95.0)
        pm3 = ProteinHit([a1, b2], 1, 6, 1, 18)
        self.assertNotEqual(pm1.protein_hit_id, pm3.protein_hit_id)

    def test_extra_match_in_middle_changes_protein_id(self):
        # Two identical ProteinHit objects should have the same ID
        a1 = Match("QID","TID",1,3,1,9,0.0,100.0)
        a2 = Match("QID","TID",4,6,10,18,1e-12,95.0)
        pm1 = ProteinHit([a1, a2], 1, 6, 1, 18)
        pm2 = ProteinHit([a1, a2], 1, 6, 1, 18)
        self.assertEqual(pm1.protein_hit_id, pm2.protein_hit_id)

        # Add new Match in middle, changes protein id
        b2 = Match("QID","TID",3,4,9,10,1e-12,95.0)
        pm3 = ProteinHit([a1, b2, a2], 1, 6, 1, 18)
        self.assertNotEqual(pm1.protein_hit_id, pm3.protein_hit_id)

    def test_target_sequence_and_target_sequence_translated_on_fwd_strand(self):

        # fwd
        m1 = self.makeM(34, 48, 1001, 1020, "Q", "S1", target_sequence = "ATGATGATG")
        self.assertEqual(m1.target_sequence_translated(), "MMM")

        # rev, but use target_sequence which is already reversed
        m2 = self.makeM(25, 38, 1201, 1120, "Q", "S1", target_sequence = "ATGATGATG")
        self.assertEqual(m2.target_sequence_translated(), "MMM")

    def test_translation_and_collated_protein_sequence(self):
        m1 = self.makeM(1, 3, 11, 19, "Q", "S1", target_sequence = "ATGGAATTT")
        m2 = self.makeM(2, 4, 22, 30, "Q", "S1", target_sequence = "GAAGTGGGG")

        self.assertEqual(m1.target_sequence_translated(), "MEF")
        self.assertEqual(m2.target_sequence_translated(), "EVG")

        pm = ProteinHit([m1, m2], 1, 4, 11, 30)
        collated = pm.collated_protein_sequence
        self.assertEqual(collated, "M(EF/EV)G")

    def test_collated_protein_sequence_does_not_include_leading_Xs(self):
        m1 = self.makeM(2, 4, 11, 19, "Q", "S1", target_sequence = "ATGGAATTT")
        m2 = self.makeM(5, 7, 22, 30, "Q", "S1", target_sequence = "GAAGTGGGG")

        self.assertEqual(m1.target_sequence_translated(), "MEF")
        self.assertEqual(m2.target_sequence_translated(), "EVG")

        pm = ProteinHit([m1, m2], 2, 7, 11, 30)
        collated = pm.collated_protein_sequence
        self.assertEqual(collated, "MEFEVG")

    def test_collate_handles_single_match(self):
        a = Match("Q","T",1,3,1,9,0.0,100.0,target_sequence="ATGGAATTT")    # MEF
        pm = ProteinHit([a],1,3,1,9)
        collated = pm.collated_protein_sequence
        self.assertEqual(collated, "MEF")

    def test_collate_handles_gaps_and_overlaps(self):
        a = Match("Q","T",1,3,1,9,0.0,100.0,target_sequence="ATGGAATTT")    # MEF
        b = Match("Q","T",3,5,10,18,0.0,100.0,target_sequence="GAAGTGGGG")  # EVG
        c = Match("Q","T",9,9,30,32,0.0,100.0,target_sequence="ATG")        # M
        pm = ProteinHit([a,b,c],1,9,1,32)
        collated = pm.collated_protein_sequence
        self.assertEqual(collated, "ME(F/E)VGXXXM")

    def test_collate_handles_gaps_within_match(self):
        a = Match("Q","T",1,3,1,6,0.0,100.0,target_sequence="ATGGAA")       # ME - but matching to 3 bps of query
        b = Match("Q","T",3,5,10,18,0.0,100.0,target_sequence="GAAGTGGGG")  # EVG
        c = Match("Q","T",9,9,30,32,0.0,100.0,target_sequence="ATG")        # M
        pm = ProteinHit([a,b,c],1,9,1,32)
        collated = pm.collated_protein_sequence
        self.assertEqual(collated, "M(E/E)VGXXXM")

    def test_collate_handles_insertions_within_match(self):
        a = Match("Q","T",1,3,1,12,0.0,100.0,target_sequence="ATGGAATTTTTT") # MEFF - but matching to 3 bps of query
        b = Match("Q","T",3,5,10,18,0.0,100.0,target_sequence="GAAGTGGGG")   # EVG
        c = Match("Q","T",9,9,30,32,0.0,100.0,target_sequence="ATG")         # M
        pm = ProteinHit([a,b,c],1,9,1,32)
        collated = pm.collated_protein_sequence
        self.assertEqual(collated, "MEF(F/E)VGXXXM")
