import unittest
from needle.expand import find_more_matches_at_locus


class TestLookingForProtein(unittest.TestCase):

    def make_hmm_row(self, evalue, score, hmm_from, hmm_to, ali_from, ali_to, matched_sequence):
        return dict( 
            target_name = "target",
            target_accession = "target",
            query_name = "query",
            query_accession = "query",
            seq_evalue = evalue,
            seq_score = score,
            dom_evalue = evalue,
            dom_evalue_cond = evalue,
            dom_score = score,
            query_length = abs(hmm_from-hmm_to)+1,
            hmm_from = hmm_from,
            hmm_to = hmm_to,
            target_length = abs(ali_from-ali_to)+1,
            ali_from = ali_from,
            ali_to = ali_to,
            matched_sequence = matched_sequence
        )

    def test_finds_original_protein_on_fwd_strand_if_there_are_no_more_matches(self):

        rows = [
            self.make_hmm_row(0.01, 10,  1, 10, 101, 130, "K"*10),
            self.make_hmm_row(0.01, 10, 11, 20, 201, 230, "K"*10),
        ]

        matches = find_more_matches_at_locus(None, None, 101, 230, "A"*500, "target", None, None, 1, hmm_rows = rows)
        self.assertEqual(len(matches), 2)
        self.assertEqual([m.target_start for m in matches], [101, 201])

    def test_finds_original_protein_on_rev_strand_if_there_are_no_more_matches(self):

        rows = [
            self.make_hmm_row(0.01, 10, 11, 20, 130, 101, "K"*10),
            self.make_hmm_row(0.01, 10,  1, 10, 230, 201, "K"*10),
        ]

        matches = find_more_matches_at_locus(None, None, 230, 101, "T"*500, "target", None, None, -1, hmm_rows = rows)
        self.assertEqual(len(matches), 2)
        self.assertEqual([m.target_start for m in matches], [230, 130])

    def test_expand_upstream_on_fwd_strand(self):

        rows = [
            self.make_hmm_row(0.01, 10,  1, 10, 101, 130, "K"*10),
            self.make_hmm_row(0.01, 10, 11, 20, 201, 230, "K"*10),
            self.make_hmm_row(0.01, 10, 21, 30, 301, 330, "K"*10),
        ]

        matches = find_more_matches_at_locus(None, None, 201, 330, "A"*500, "target", None, None, 1, hmm_rows = rows)
        self.assertEqual(len(matches), 3)
        self.assertEqual([m.target_start for m in matches], [101, 201, 301])

    def test_expand_downstream_on_fwd_strand(self):

        rows = [
            self.make_hmm_row(0.01, 10,  1, 10, 101, 130, "K"*10),
            self.make_hmm_row(0.01, 10, 11, 20, 201, 230, "K"*10),
            self.make_hmm_row(0.01, 10, 21, 30, 301, 330, "K"*10),
        ]

        matches = find_more_matches_at_locus(None, None, 101, 230, "A"*500, "target", None, None, 1, hmm_rows = rows)
        self.assertEqual(len(matches), 3)
        self.assertEqual([m.target_start for m in matches], [101, 201, 301])

    def test_expand_upstream_on_rev_strand(self):

        rows = [
            self.make_hmm_row(0.01, 10, 21, 30, 130, 101, "K"*10),
            self.make_hmm_row(0.01, 10, 11, 20, 230, 201, "K"*10),
            self.make_hmm_row(0.01, 10,  1, 10, 330, 301, "K"*10),
        ]

        matches = find_more_matches_at_locus(None, None, 230, 101, "T"*500, "target", None, None, -1, hmm_rows = rows)
        self.assertEqual(len(matches), 3)
        self.assertEqual([m.target_start for m in matches], [330, 230, 130])

    def test_expand_downstream_on_rev_strand(self):

        rows = [
            self.make_hmm_row(0.01, 10, 21, 30, 130, 101, "K"*10),
            self.make_hmm_row(0.01, 10, 11, 20, 230, 201, "K"*10),
            self.make_hmm_row(0.01, 10,  1, 10, 330, 301, "K"*10),
        ]

        matches = find_more_matches_at_locus(None, None, 330, 201, "T"*500, "target", None, None, -1, hmm_rows = rows)
        self.assertEqual(len(matches), 3)
        self.assertEqual([m.target_start for m in matches], [330, 230, 130])

    def test_can_expand_both_directions(self):

        rows = [
            self.make_hmm_row(0.01, 10,  1, 10, 101, 130, "K"*10),
            self.make_hmm_row(0.01, 10, 11, 20, 201, 230, "K"*10),
            self.make_hmm_row(0.01, 10, 21, 30, 301, 330, "K"*10),
            self.make_hmm_row(0.01, 10, 31, 40, 401, 430, "K"*10),
        ]

        matches = find_more_matches_at_locus(None, None, 201, 330, "A"*500, "target", None, None, 1, hmm_rows = rows)
        self.assertEqual(len(matches), 4)
        self.assertEqual([m.target_start for m in matches], [101, 201, 301, 401])

    def test_does_not_expand_into_potential_overlapping_genes(self):

        rows = [
            self.make_hmm_row(0.01, 10, 11, 20, 101, 130, "K"*10),  # overlaps with next one
            self.make_hmm_row(0.01, 10, 11, 20, 201, 230, "K"*10),
            self.make_hmm_row(0.01, 10, 21, 30, 301, 330, "K"*10),
            self.make_hmm_row(0.01, 10, 31, 40, 401, 430, "K"*10),
        ]

        matches = find_more_matches_at_locus(None, None, 201, 330, "A"*500, "target", None, None, 1, hmm_rows = rows)
        self.assertEqual(len(matches), 3)
        self.assertEqual([m.target_start for m in matches], [201, 301, 401])

        rows = [
            self.make_hmm_row(0.01, 10,  1, 10, 101, 130, "K"*10),
            self.make_hmm_row(0.01, 10, 11, 20, 201, 230, "K"*10),
            self.make_hmm_row(0.01, 10, 21, 30, 301, 330, "K"*10),
            self.make_hmm_row(0.01, 10, 21, 30, 401, 430, "K"*10),  # overlaps with prev one
        ]

        matches = find_more_matches_at_locus(None, None, 201, 330, "A"*500, "target", None, None, 1, hmm_rows = rows)
        self.assertEqual(len(matches), 3)
        self.assertEqual([m.target_start for m in matches], [101, 201, 301])
