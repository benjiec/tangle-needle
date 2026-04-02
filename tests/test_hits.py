import os
import io
import tempfile
import unittest

from needle.match import order_matches_for_junctions

from needle.hits import (
    generate_transition_candidates,
    score_and_select_best_transition,
    stitch_cleaned_sequence,
    Candidate,
    hmm_clean_protein,
    hmm_clean,
    adjust_target_coordinates,
    find_matches_at_locus
)
import needle.hits as hits_mod

from needle.detect import Results
from needle.match import group_matches, ProteinHit, Match
from needle.hits import write_fasta_record, export_protein_hits


class TestCleaningSequenceWithHMM(unittest.TestCase):

    def test_generate_transition_candidates_overlap(self):
        left = "ABCDEFX"
        right = "yefghij"

        cands = generate_transition_candidates(left, right, overlap_len=4, gap_len=0, overlap_flanking_len=2)
        self.assertEqual(len(cands), 5)

        self.assertEqual(
            [c.left_trimmed for c in cands],
            [
                4,     # k=0: left trims 4, right trims 0
                3,     # k=1: left trims 3, right trims 1
                2,     # k=2: left trims 2, right trims 2
                1,     # k=3: left trims 1, right trims 3
                0,     # k=4: left trims 0, right trims 4
            ],
        )
        self.assertEqual(
            [c.assigned_overlap_to_left for c in cands],
            [
                0,     # k=0: left trims 4, right trims 0
                1,     # k=1: left trims 3, right trims 1
                2,     # k=2: left trims 2, right trims 2
                3,     # k=3: left trims 1, right trims 3
                4,     # k=4: left trims 0, right trims 4
            ],
        )
        self.assertEqual(
            [c.stitched for c in cands],
            [
                "ABCyefghij",     # k=0: left trims 4, right trims 0
                "ABCDefghij",     # k=1: left trims 3, right trims 1
                "ABCDEfghij",     # k=2: left trims 2, right trims 2
                "ABCDEFghij",     # k=3: left trims 1, right trims 3
                "ABCDEFXhij",     # k=4: left trims 0, right trims 4
            ],
        )
        self.assertEqual(
            [c.right_kept for c in cands],
            [
                "yefghij",     # k=0: left trims 4, right trims 0
                "efghij",     # k=1: left trims 3, right trims 1
                "fghij",     # k=2: left trims 2, right trims 2
                "ghij",     # k=3: left trims 1, right trims 3
                "hij",     # k=4: left trims 0, right trims 4
            ],
        )
        self.assertEqual(
            [c.window_seq for c in cands],
            [
                "BCyefghi",     # k=0: left trims 4, right trims 0
                "BCDefghi",     # k=1: left trims 3, right trims 1
                "BCDEfghi",     # k=2: left trims 2, right trims 2
                "BCDEFghi",     # k=3: left trims 1, right trims 3
                "BCDEFXhi",     # k=4: left trims 0, right trims 4
            ],
        )

    def test_generate_transition_candidates_overlap_same_as_full(self):
        left = "DEFX"
        right = "yefg"

        cands = generate_transition_candidates(left, right, overlap_len=4, gap_len=0, overlap_flanking_len=2)
        self.assertEqual(len(cands), 5)

        self.assertEqual(
            [c.left_trimmed for c in cands],
            [
                4,     # k=0: left trims 4, right trims 0
                3,     # k=1: left trims 3, right trims 1
                2,     # k=2: left trims 2, right trims 2
                1,     # k=3: left trims 1, right trims 3
                0,     # k=4: left trims 0, right trims 4
            ],
        )
        self.assertEqual(
            [c.stitched for c in cands],
            [
                "yefg",     # k=0: left trims 4, right trims 0
                "Defg",     # k=1: left trims 3, right trims 1
                "DEfg",     # k=2: left trims 2, right trims 2
                "DEFg",     # k=3: left trims 1, right trims 3
                "DEFX",     # k=4: left trims 0, right trims 4
            ],
        )
        self.assertEqual(
            [c.right_kept for c in cands],
            [
                "yefg",     # k=0: left trims 4, right trims 0
                "efg",     # k=1: left trims 3, right trims 1
                "fg",     # k=2: left trims 2, right trims 2
                "g",     # k=3: left trims 1, right trims 3
                "",     # k=4: left trims 0, right trims 4
            ],
        )
        self.assertEqual(
            [c.window_seq for c in cands],
            [
                "yefg",     # k=0: left trims 4, right trims 0
                "Defg",     # k=1: left trims 3, right trims 1
                "DEfg",     # k=2: left trims 2, right trims 2
                "DEFg",     # k=3: left trims 1, right trims 3
                "DEFX",     # k=4: left trims 0, right trims 4
            ],
        )
  
    def test_generate_transition_candidates_gap(self):
        cands_gap = generate_transition_candidates("AAA", "bbb", overlap_len=0, gap_len=2, overlap_flanking_len=2)
        self.assertEqual(len(cands_gap), 1)
        self.assertEqual(cands_gap[0].assigned_overlap_to_left, None)
        self.assertEqual(cands_gap[0].left_trimmed, 0)
        self.assertEqual(cands_gap[0].stitched, "AAAXXbbb")
        self.assertEqual(cands_gap[0].right_kept, "XXbbb")
        self.assertEqual(cands_gap[0].window_seq, "AAXXbb")

    @staticmethod
    def makeM(query_start, query_end, target_start, target_end):
        return Match(
            query_accession=None, target_accession=None, e_value=0, identity=None,
            query_start=query_start, query_end=query_end, target_start=target_start, target_end=target_end,
            target_sequence="A"*(abs(target_start-target_end)+1))

    def test_stitch_cleaned_sequence_basic(self):
        left = self.makeM(1, 5, 1, 15)
        right = self.makeM(4, 8, 10, 24)
        aa_map = {id(left):"ABCDE", id(right):"DEFGH"}

        pairs = order_matches_for_junctions([left, right])  # type: ignore
        self.assertEqual(pairs[0][2], 2)
        cand = Candidate(assigned_overlap_to_left=1, window_seq="", stitched="ABCDEFGH", left_trimmed=2, right_kept="DEFGH")
        stitched = stitch_cleaned_sequence([(left, right, None, None)], {0: cand}, aa_map)  # type: ignore
        self.assertEqual(stitched, "ABCDEFGH")

    def test_stitch_cleaned_sequence_multiple_blocks_mixed(self):
        a = self.makeM(1, 5, 1, 15)
        b = self.makeM(4, 9, 10, 27)
        c = self.makeM(13, 15, 37, 45)
        aa_map = {id(a):"ABCDE", id(b):"DEFGHI", id(c):"KLM"}

        pairs = order_matches_for_junctions([a,b,c])  # type: ignore
        cand0 = Candidate(assigned_overlap_to_left=1, window_seq="", stitched="ABCDEFGHI", left_trimmed=1, right_kept="EFGHI")
        cand1 = Candidate(assigned_overlap_to_left=None, window_seq="", stitched="DEFGHIXXXKLM", left_trimmed=0, right_kept="XXXKLM")
        stitched = stitch_cleaned_sequence(
          [(a,b, None, None), (b,c, None, None)],
          {0:cand0, 1:cand1}, aa_map)
        self.assertEqual(stitched, "ABCDEFGHIXXXKLM")

    def test_hmm_cleaned_protein_integration_with_mock_scoring(self):
        a = Match("Q","T",1,3,1,9,0.0,100.0,False); a.target_sequence="ATGGAATTT"    # MEF
        b = Match("Q","T",3,5,10,18,0.0,100.0,False); b.target_sequence="GAAGTGGGG"  # EVG
        c = Match("Q","T",9,9,30,32,0.0,100.0,False); c.target_sequence="ATG"        # M
        pm = ProteinHit([a,b,c],1,9,1,32)

        orig = hits_mod.score_and_select_best_transition
        def _fake(cands, hmm): 
            for x in cands:
                if x.assigned_overlap_to_left == 1: return x
            return cands[0]
        try:
            hits_mod.score_and_select_best_transition = _fake
            cleaned_pm = hmm_clean_protein(pm, "dummy.hmm", overlap_flanking_len=5, min_query_match_len=0)
            cleaned = cleaned_pm.collated_protein_sequence
        finally:
            hits_mod.score_and_select_best_transition = orig

        self.assertEqual(cleaned, "MEFVGXXXM")

    def test_hmm_clean_protein_adjusts_overlap_coordinates(self):
        # Overlap: a(1..5), b(4..9) => overlap 2; choose k=1; c(13..15) should shift by 1

        a = Match("Q","T",1,5,1,15,0.0,100.0,False); a.target_sequence="ATG"*5         # 'M'*5
        b = Match("Q","T",4,9,16,33,0.0,100.0,False); b.target_sequence="GAA"*6        # 'E'*6
        c = Match("Q","T",13,17,40,51,0.0,100.0,False); c.target_sequence="ATG"*5      # 'M'*5
        pm = ProteinHit([a,b,c],1,15,1,51)

        orig = hits_mod.score_and_select_best_transition
        def _fake(cands, hmm):
            for x in cands:
                if x.assigned_overlap_to_left == 1:
                    return x
            return cands[0]
        try:
            hits_mod.score_and_select_best_transition = _fake  # type: ignore
            cleaned_pm = hmm_clean_protein(pm, "dummy.hmm", overlap_flanking_len=5)
        finally:
            hits_mod.score_and_select_best_transition = orig  # type: ignore

        nm = cleaned_pm.matches
        self.assertEqual(len(nm), 3)
        # a.end reduced by 1
        self.assertEqual(nm[0].query_start, 1)
        self.assertEqual(nm[0].query_end, 4)
        # b.start becomes a.end+1 = 5; b end remains 9 (no shift of downstream blocks)
        self.assertEqual(nm[1].query_start, 5)
        self.assertEqual(nm[1].query_end, 9)
        # c unchanged (no shifting of downstream matches)
        self.assertEqual(nm[2].query_start, 13)
        self.assertEqual(nm[2].query_end, 17)

    def test_adjust_target_coordinates_gap_keeps_blocks(self):
        # query acc, target acc, query start, query end, target start, target end
        left = Match("q","t",1,5,100,114,0.0,100.0,False); left.target_sequence="A"*15
        right = Match("q","t",8,12,200,214,0.0,100.0,False); right.target_sequence="C"*15
        cand = Candidate(assigned_overlap_to_left=None, window_seq="", stitched="", left_trimmed=None, right_kept="")
        new_left, new_right = adjust_target_coordinates(left, right, cand)
        self.assertEqual((new_left.query_start, new_left.query_end), (1,5))
        self.assertEqual((new_right.query_start, new_right.query_end), (8,12))
        self.assertEqual(new_left.target_sequence, "A"*15)
        self.assertEqual(new_right.target_sequence, "C"*15)

    def test_adjust_target_coordinates_overlap_k0_no_change(self):
        # query acc, target acc, query start, query end, target start, target end
        left = Match("q","t",1,5,100,114,0.0,100.0,False); left.target_sequence="A"*15
        right = Match("q","t",5,9,200,214,0.0,100.0,False); right.target_sequence="C"*15
        cand = Candidate(assigned_overlap_to_left=0, window_seq="", stitched="", left_trimmed=1, right_kept="")
        new_left, new_right = adjust_target_coordinates(left, right, cand)
        self.assertEqual((new_left.query_start, new_left.query_end), (1,4))
        self.assertEqual((new_right.query_start, new_right.query_end), (5,9))
        self.assertEqual(new_left.target_sequence, "A"*12)
        self.assertEqual(new_right.target_sequence, "C"*15)

    def test_adjust_target_coordinates_overlap_k1_trims_and_adjacent(self):
        # query acc, target acc, query start, query end, target start, target end
        left = Match("q","t",1,5,100,114,0.0,100.0,False); left.target_sequence="A"*15
        right = Match("q","t",4,9,200,217,0.0,100.0,False); right.target_sequence="C"*18
        cand = Candidate(assigned_overlap_to_left=1, window_seq="", stitched="", left_trimmed=1, right_kept="")
        new_left, new_right = adjust_target_coordinates(left, right, cand)
        self.assertEqual((new_left.query_start, new_left.query_end), (1,4))
        self.assertEqual((new_right.query_start, new_right.query_end), (5,9))
        self.assertEqual(new_left.target_sequence, "A"*12)
        self.assertEqual(new_right.target_sequence, "C"*15)

    def test_adjust_target_coordinates_overlap_trim_all_from_left(self):
        # query acc, target acc, query start, query end, target start, target end
        left = Match("q","t",1,1,100,102,0.0,100.0,False); left.target_sequence="A"*3
        right = Match("q","t",1,4,200,211,0.0,100.0,False); right.target_sequence="C"*12
        cand = Candidate(assigned_overlap_to_left=0, window_seq="", stitched="", left_trimmed=1, right_kept="")
        new_left, new_right = adjust_target_coordinates(left, right, cand)
        self.assertEqual((new_left.query_start, new_left.query_end), (1,0))
        self.assertEqual((new_right.query_start, new_right.query_end), (1,4))
        self.assertEqual(new_left.target_sequence, "")
        self.assertEqual(new_right.target_sequence, "C"*12)


class TestRefiningHitsWithHMM(unittest.TestCase):

    def test_find_matches_at_locus_incrementally_search_for_more_matches_to_specified_range_beyond_last_position(self):

        orig = hits_mod.hmm_search_genome

        searched = []

        def fake_hmm_search_genome(_hmm_file, _gs, target_accession, target_left, target_right, strand):
            first_match = [
                dict(target_name="cand_0", dom_score=100, dom_evalue=0.0001, hmm_from=5, hmm_to=10, ali_from=10001, ali_to=10018, matched_sequence="F"*6),
                dict(target_name="cand_1", dom_score=100, dom_evalue=0.0002, hmm_from=9, hmm_to=15, ali_from=11021, ali_to=11041, matched_sequence="F"*7),
                dict(target_name="cand_2", dom_score=100, dom_evalue=0.0003, hmm_from=16, hmm_to=20, ali_from=11051, ali_to=11065, matched_sequence="F"*5)
            ]
            second_match = [
                # add a new match
                dict(target_name="cand_0", dom_score=100, dom_evalue=0.0001, hmm_from=1, hmm_to=4, ali_from=9001, ali_to=9012, matched_sequence="F"*4),
                dict(target_name="cand_0", dom_score=100, dom_evalue=0.0001, hmm_from=5, hmm_to=10, ali_from=10001, ali_to=10018, matched_sequence="F"*6),
                dict(target_name="cand_1", dom_score=100, dom_evalue=0.0002, hmm_from=9, hmm_to=15, ali_from=11021, ali_to=11041, matched_sequence="F"*7),
                dict(target_name="cand_2", dom_score=100, dom_evalue=0.0003, hmm_from=16, hmm_to=20, ali_from=11051, ali_to=11065, matched_sequence="F"*5)
            ]

            searched.append((target_left, target_right))

            if target_left == 10001: # initial
                return first_match
            elif target_left < 10001:
                return second_match

        try:
            hits_mod.hmm_search_genome = fake_hmm_search_genome
          
            old_matches = [
                Match(query_accession="Q", target_accession="T", query_start=5, query_end=10, target_start=10001, target_end=10018, e_value=0.1, identity=None)
            ]

            new_matches = find_matches_at_locus(
                old_matches,
                "T"*15000,
                10001, 12000, "hmmfile", step=2000, max_search_distance=6000
            )

            self.assertNotEqual(new_matches, None)
            self.assertEqual(len(new_matches), 4)

            # last found match is at 9001, so we can search to 3001
            self.assertEqual(searched, [(10001, 12000), (8001, 14000), (6001, 15000), (4001, 15000), (2001, 15000)])

            self.assertEqual(new_matches[0].query_accession, "Q")
            self.assertEqual(new_matches[0].target_accession, "T")

            self.assertEqual(new_matches[0].query_start, 1)
            self.assertEqual(new_matches[0].query_end, 4)
            self.assertEqual(new_matches[0].target_start, 9001)
            self.assertEqual(new_matches[0].target_end, 9012)

            self.assertEqual(new_matches[1].query_start, 5)
            self.assertEqual(new_matches[1].query_end, 10)
            self.assertEqual(new_matches[1].target_start, 10001)
            self.assertEqual(new_matches[1].target_end, 10018)

            self.assertEqual(new_matches[2].query_start, 9)
            self.assertEqual(new_matches[2].query_end, 15)
            self.assertEqual(new_matches[2].target_start, 11021)
            self.assertEqual(new_matches[2].target_end, 11041)

            self.assertEqual(new_matches[3].query_start, 16)
            self.assertEqual(new_matches[3].query_end, 20)
            self.assertEqual(new_matches[3].target_start, 11051)
            self.assertEqual(new_matches[3].target_end, 11065)

        finally:
            hits_mod.hmm_search_genome = orig

    def test_find_matches_at_locus_incrementally_search_on_rev_strand_as_well(self):

        orig = hits_mod.hmm_search_genome

        searched = []
        def fake_hmm_search_genome(_hmm_file, _gs, target_accession, target_left, target_right, strand):
            first_match = [
                dict(target_name="cand_0", dom_score=100, dom_evalue=0.0001, hmm_from=5, hmm_to=10, ali_from=11018, ali_to=11001, matched_sequence="F"*6),
                dict(target_name="cand_1", dom_score=100, dom_evalue=0.0002, hmm_from=9, hmm_to=15, ali_from=10841, ali_to=10821, matched_sequence="F"*7),
                dict(target_name="cand_2", dom_score=100, dom_evalue=0.0003, hmm_from=16, hmm_to=20, ali_from=10065, ali_to=10051, matched_sequence="F"*5)
            ]

            searched.append((target_left, target_right))

            if target_right == 12000: # initial
                return first_match
            elif target_right > 12000:
                return first_match

        try:
            hits_mod.hmm_search_genome = fake_hmm_search_genome
          
            old_matches = [
                Match(query_accession="Q", target_accession="T", query_start=5, query_end=10, target_start=10018, target_end=10001, e_value=0.1, identity=None)
            ]

            new_matches = find_matches_at_locus(
                old_matches,
                "A"*25000,
                12000, 10001, "hmmfile", step=2000, max_search_distance=6000
            )

            self.assertNotEqual(new_matches, None)
            self.assertEqual(len(new_matches), 3)
            # last match is 11018, so can go to 17018
            self.assertEqual(searched, [(10001, 12000), (8001, 14000), (6001, 16000), (4001, 18000)])

        finally:
            hits_mod.hmm_search_genome = orig

    def test_find_matches_at_locus_ensures_hmm_matched_sequence_matches_translated_sequence(self):

        orig = hits_mod.hmm_search_genome
        expected_aa = "F"*6

        def fake_hmm_search_genome(_hmm_file, _gs, target_accession, target_left, target_right, strand):
            first_match = [
                dict(target_name="cand_0", dom_score=100, dom_evalue=0.1, hmm_from=5, hmm_to=10, ali_from=10001, ali_to=10018, matched_sequence=expected_aa)
            ]
            return first_match

        try:
            hits_mod.hmm_search_genome = fake_hmm_search_genome
          
            old_matches = [
                Match(query_accession="Q", target_accession="T", query_start=5, query_end=10, target_start=10001, target_end=10018, e_value=0.1, identity=None)
            ]

            # target sequence on DNA matches what HMM says
            new_matches = find_matches_at_locus(
                old_matches,
                "C"*10000+"T"*18+"C"*20000,
                10001, 20000, "hmmfile", step=2000, force_extend=True
            )
            self.assertEqual(len(new_matches), 1)
            self.assertEqual(new_matches[0].target_sequence_translated(), expected_aa)

            # target sequence does not match what HMM says, for some odd reason
            with self.assertRaises(AssertionError):
                new_matches = find_matches_at_locus(
                    old_matches,
                    # not 18 T's
                    "C"*10000+"T"*14+"C"*20000,
                    10001, 20000, "hmmfile", step=2000, force_extend=True
                )

        finally:
            hits_mod.hmm_search_genome = orig

    def test_find_matches_at_locus_stops_searching_at_boundaries(self):

        orig = hits_mod.hmm_search_genome

        searched = []
        def fake_hmm_search_genome(_hmm_file, _gs, target_accession, target_left, target_right, strand):
            first_match = [
                dict(target_name="cand_0", dom_score=100, dom_evalue=0.0001, hmm_from=5, hmm_to=10, ali_from=10001, ali_to=10018, matched_sequence="F"*6),
                dict(target_name="cand_1", dom_score=100, dom_evalue=0.0002, hmm_from=9, hmm_to=15, ali_from=11021, ali_to=11041, matched_sequence="F"*7),
                dict(target_name="cand_2", dom_score=100, dom_evalue=0.0003, hmm_from=16, hmm_to=20, ali_from=11051, ali_to=11065, matched_sequence="F"*5)
            ]
            second_match = [
                # add a new match
                dict(target_name="cand_0", dom_score=100, dom_evalue=0.0001, hmm_from=1, hmm_to=4, ali_from=9001, ali_to=9012, matched_sequence="F"*4),
                dict(target_name="cand_0", dom_score=100, dom_evalue=0.0001, hmm_from=5, hmm_to=10, ali_from=10001, ali_to=10018, matched_sequence="F"*6),
                dict(target_name="cand_1", dom_score=100, dom_evalue=0.0002, hmm_from=9, hmm_to=15, ali_from=11021, ali_to=11041, matched_sequence="F"*7),
                dict(target_name="cand_2", dom_score=100, dom_evalue=0.0003, hmm_from=16, hmm_to=20, ali_from=11051, ali_to=11065, matched_sequence="F"*5)
            ]

            searched.append((target_left, target_right))

            if target_left == 10001: # initial
                return first_match
            elif target_left < 10001:
                return second_match

        try:
            hits_mod.hmm_search_genome = fake_hmm_search_genome
          
            old_matches = [
                Match(query_accession="Q", target_accession="T", query_start=5, query_end=10, target_start=10001, target_end=10018, e_value=0.1, identity=None)
            ]

            new_matches = find_matches_at_locus(
                old_matches,
                "T"*25000,
                10001, 12000, "hmmfile", step=6000, max_search_distance=15000
            )

            self.assertEqual(searched, [(10001, 12000), (4001, 18000), (1, 24000), (1, 25000)])

        finally:
            hits_mod.hmm_search_genome = orig


class TestIO(unittest.TestCase):

    def test_write_fasta_record(self):
        a = Match("QX","TX",1,3,1,9,0.0,100.0,False); a.target_sequence="ATGGAATTT"
        pm = ProteinHit([a],1,3,1,9)
        pid = pm.protein_hit_id
        buf = io.StringIO()
        write_fasta_record(buf, pm)
        s = buf.getvalue().splitlines()
        self.assertEqual(s[0], f">{pid}")
        self.assertTrue(s[1].startswith("MEF"))

    def test_export_protein_hits_filters_and_writes(self):
        # pm1: single block => eligible
        a = Match("Q1","T1",1,3,1,9,0.0,100.0,False); a.target_sequence="ATGGAATTT"
        pm1 = ProteinHit([a],1,3,1,9, hmm_file="ignored.hmm")
        # pm2: overlapping blocks => not single sequence
        b1 = Match("Q2","T2",1,3,1,9,0.0,100.0,False); b1.target_sequence="ATGGAATTT"
        b2 = Match("Q2","T2",3,5,10,18,0.0,100.0,False); b2.target_sequence="GAAGTGGGG"
        pm2 = ProteinHit([b1,b2],1,5,1,18, hmm_file="ignored.hmm")
        with tempfile.TemporaryDirectory() as d:
            p1 = os.path.join(d, "prot.faa")
            p2 = os.path.join(d, "prot.tsv")
            export_protein_hits("GENOMEZ", [pm1, pm2], p1, p2)
            with open(p1) as f:
                lines = [l.strip() for l in f if l.strip()]
            with open(p2) as f:
                lines2 = [l.strip() for l in f if l.strip()]
            # fasta: 2 lines (header+seq)
            self.assertEqual(len(lines), 2)
            # nuc rows: only pm1 yields rows => 1 line + header
            self.assertEqual(len(lines2), 2)

    def test_export_protein_hits_appends_not_overwrites(self):
        a = Match("QX","TX",1,3,1,9,0.0,100.0,False); a.target_sequence="ATGGAATTT"
        pm = ProteinHit([a],1,3,1,9, hmm_file="ignored.hmm")
        with tempfile.TemporaryDirectory() as d:
            p1 = os.path.join(d, "prot.faa")
            p2 = os.path.join(d, "nuc.tsv")
            export_protein_hits("GENOME1", [pm], p1, p2)
            # capture initial counts
            with open(p1) as f:
                lines_prot_1 = [l for l in f if l.strip()]
            with open(p2) as f:
                lines_nuc_1 = [l for l in f if l.strip()]
            # second export appends
            export_protein_hits("GENOME1", [pm], p1, p2)
            with open(p1) as f:
                lines_prot_2 = [l for l in f if l.strip()]
            with open(p2) as f:
                lines_nuc_2 = [l for l in f if l.strip()]
            # proteins.tsv: header + 1 row initially; after append expect +2 lines
            self.assertEqual(len(lines_prot_1), 2)
            self.assertEqual(len(lines_prot_2), 4)
            # nuc.tsv: header + 1 row initially; after append expect +1 row
            self.assertEqual(len(lines_nuc_1), 2)
            self.assertEqual(len(lines_nuc_2), 3)
