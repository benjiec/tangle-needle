import unittest
from needle.expand import hmm_expand_protein
from needle.match import ProteinHit, Match


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

        a = Match("Q", "T", 11, 20, 101, 130, 0.0, 100.0)
        b = Match("Q", "T", 21, 30, 201, 230, 0.0, 100.0)

        rows = [
            self.make_hmm_row(0.01, 10, 11, 20, 101, 130, "K"*10),
            self.make_hmm_row(0.01, 10, 21, 30, 201, 230, "K"*10),
        ]

        proteins = hmm_expand_protein([a, b], {"T": "A"*500}, None, threshold = None, hmm_rows = rows)
        self.assertEqual(len(proteins), 1)
        self.assertEqual([m.target_start for m in proteins[0].matches], [101, 201])

    def test_finds_expanded_protein_on_fwd_strand(self):

        a = Match("Q", "T", 11, 20, 101, 130, 0.0, 100.0)
        b = Match("Q", "T", 21, 30, 201, 230, 0.0, 100.0)

        rows = [
            self.make_hmm_row(0.01, 10, 11, 20, 101, 130, "K"*10),
            self.make_hmm_row(0.01, 10,  1, 10,  11,  40, "K"*10),
            self.make_hmm_row(0.01, 10, 21, 30, 201, 230, "K"*10),
            self.make_hmm_row(0.01, 10, 31, 40, 401, 430, "K"*10),
        ]

        proteins = hmm_expand_protein([a, b], {"T": "A"*500}, None, threshold = None, hmm_rows = rows)

        self.assertEqual(len(proteins), 1)
        self.assertEqual(proteins[0].on_reverse_strand, False)
        self.assertEqual(proteins[0].target_start, 11)
        self.assertEqual(proteins[0].target_end, 430)
        self.assertEqual(len(proteins[0].matches), 4)

    def test_finds_expanded_protein_on_rev_strand(self):
        a = Match("Q", "T", 11, 20, 230, 201, 0.0, 100.0)
        b = Match("Q", "T", 21, 30, 130, 101, 0.0, 100.0)

        rows = [
            self.make_hmm_row(0.01, 10, 21, 30, 130, 101, "F"*10),
            self.make_hmm_row(0.01, 10, 31, 40,  40,  11, "F"*10),
            self.make_hmm_row(0.01, 10, 11, 20, 230, 201, "F"*10),
            self.make_hmm_row(0.01, 10,  1, 10, 430, 401, "F"*10),
        ]

        proteins = hmm_expand_protein([a, b], {"T": "A"*500}, None, threshold = None, hmm_rows = rows)
        self.assertEqual(len(proteins), 1)

        self.assertEqual(proteins[0].on_reverse_strand, True)
        self.assertEqual(proteins[0].target_start, 430)
        self.assertEqual(proteins[0].target_end, 11)
        self.assertEqual(len(proteins[0].matches), 4)

    def test_finds_multiple_expanded_proteins_on_both_strands(self):
        a = Match("Q", "T", 11, 20, 230, 201, 0.0, 100.0)
        b = Match("Q", "T", 21, 30, 130, 101, 0.0, 100.0)

        rows = [
            self.make_hmm_row(0.01, 10, 21, 30, 130, 101, "F"*10),
            self.make_hmm_row(0.01, 10, 4,  13,  41,  70, "K"*10),
            self.make_hmm_row(0.01, 10, 31, 40,  40,  11, "F"*10),
            self.make_hmm_row(0.01, 10, 11, 20, 230, 201, "F"*10),
            self.make_hmm_row(0.01, 10,  1, 10, 430, 401, "F"*10),
        ]

        proteins = hmm_expand_protein([a, b], {"T": "A"*500}, None, threshold = None, hmm_rows = rows)
        self.assertEqual(len(proteins), 2)

        self.assertEqual(proteins[0].on_reverse_strand, True)
        self.assertEqual(proteins[0].target_start, 430)
        self.assertEqual(proteins[0].target_end, 11)
        self.assertEqual(len(proteins[0].matches), 4)

        self.assertEqual(proteins[1].on_reverse_strand, False)
        self.assertEqual(proteins[1].target_start, 41)
        self.assertEqual(proteins[1].target_end, 70)
        self.assertEqual(len(proteins[1].matches), 1)


    def test_expand_protein_rejects_protein_below_threshold(self):
        a = Match("Q", "T", 11, 20, 201, 230, 0.0, 100.0)
        b = Match("Q", "T", 21, 30, 301, 330, 0.0, 100.0)

        rows = [
            self.make_hmm_row(0.01, 10,  1, 10, 101, 130, "K"*10),
            self.make_hmm_row(0.01, 10, 11, 20, 201, 230, "K"*10),
            self.make_hmm_row(0.01, 10, 21, 30, 301, 330, "K"*10),
            self.make_hmm_row(0.01, 10, 31, 40, 401, 430, "K"*10),
        ]

        proteins = hmm_expand_protein([a, b], {"T": "A"*500}, None, threshold = None, hmm_rows = rows)
        self.assertEqual(len(proteins), 1)

        proteins = hmm_expand_protein([a, b], {"T": "A"*500}, None, threshold = 41, hmm_rows = rows)
        self.assertEqual(len(proteins), 0)

        proteins = hmm_expand_protein([a, b], {"T": "A"*500}, None, threshold = 40, hmm_rows = rows)
        self.assertEqual(len(proteins), 1)
