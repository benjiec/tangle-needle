import os
import tempfile
import unittest

from needle.hmm import parse_hmmsearch_domtbl, removable_gap_starts_here, Aligned


class TestHMMSearch(unittest.TestCase):

    def test_parse_hmmsearch_domtbl_returns_matches(self):

        with tempfile.TemporaryDirectory() as temp_dir:

            temp_file_path = os.path.join(temp_dir, 'domtbl.txt')
            with open(temp_file_path, 'w') as f:
                f.write("# target name  accession  tlen  query name  accession  qlen  E-value  score  bias  #  of  c-Evalue  i-Evalue  score  bias  from  to  from  to\n")
                f.write("cand_0 _ 81 _ _ 5 0.1 100 _ _ _ 0.005 0.05 201 _ 11 12 13 14\n")
                f.write("cand_0 _ 82 _ _ 6 0.2 101 _ _ _ 0.004 0.04 202 _ 21 22 23 24\n")
                f.close()

            matches = parse_hmmsearch_domtbl(temp_file_path)
            self.assertEqual(len(matches), 2)
            self.assertEqual(matches[0]["target_name"], "cand_0")
            self.assertEqual(matches[0]["target_length"], 81)
            self.assertEqual(matches[0]["seq_evalue"], 0.1)
            self.assertEqual(matches[0]["seq_score"], 100)
            self.assertEqual(matches[0]["dom_evalue"], 0.05) # using i-Evalue, not c-Evalue
            self.assertEqual(matches[0]["dom_evalue_cond"], 0.005)
            self.assertEqual(matches[0]["dom_score"], 201)
            self.assertEqual(matches[0]["query_length"], 5)
            self.assertEqual(matches[0]["hmm_from"], 11)
            self.assertEqual(matches[0]["hmm_to"], 12)
            self.assertEqual(matches[0]["ali_from"], 13)
            self.assertEqual(matches[0]["ali_to"], 14)
            self.assertEqual(matches[1]["target_name"], "cand_0")
            self.assertEqual(matches[1]["target_length"], 82)
            self.assertEqual(matches[1]["seq_evalue"], 0.2)
            self.assertEqual(matches[1]["seq_score"], 101)
            self.assertEqual(matches[1]["dom_evalue"], 0.04) # using i-Evalue, not c-Evalue
            self.assertEqual(matches[1]["dom_evalue_cond"], 0.004)
            self.assertEqual(matches[1]["dom_score"], 202)
            self.assertEqual(matches[1]["query_length"], 6)
            self.assertEqual(matches[1]["hmm_from"], 21)
            self.assertEqual(matches[1]["hmm_to"], 22)
            self.assertEqual(matches[1]["ali_from"], 23)
            self.assertEqual(matches[1]["ali_to"], 24)

    def test_parse_hmmsearch_domtbl_asserts_has_expected_headers(self):

        with tempfile.TemporaryDirectory() as temp_dir:

            temp_file_path = os.path.join(temp_dir, 'domtbl.txt')
            with open(temp_file_path, 'w') as f:
                # BAD unexpected header
                f.write("# target name  accession  tlen  query name  accession  qlen  E-value  score  from  to  from  to\n")
                f.write("cand_0 _ _ _ _ _ 0.1 100 11 12 13 14\n")
                f.write("cand_0 _ _ _ _ _ 0.2 101 21 22 23 24\n")
                f.close()

            with self.assertRaises(AssertionError):
                parse_hmmsearch_domtbl(temp_file_path)


class TestAlignmentParsing(unittest.TestCase):

    def test_removable_gap_starts_here_detects_gaps_beyond_tolerance(self):
        query_sequence = "fff........nnnnnn"
        hit_sequence   = "nnnnnnnnnnnnnnnnn"

        self.assertEqual(removable_gap_starts_here(query_sequence, hit_sequence, 0, gap_tolerated=7, always_remove_gap_with_star=True),
                         None)
        self.assertEqual(removable_gap_starts_here(query_sequence, hit_sequence, 1, gap_tolerated=7, always_remove_gap_with_star=True),
                         None)
        self.assertEqual(removable_gap_starts_here(query_sequence, hit_sequence, 2, gap_tolerated=7, always_remove_gap_with_star=True),
                         None)
        self.assertEqual(removable_gap_starts_here(query_sequence, hit_sequence, 3, gap_tolerated=7, always_remove_gap_with_star=True),
                         11)
        self.assertEqual(removable_gap_starts_here(query_sequence, hit_sequence, 11, gap_tolerated=7, always_remove_gap_with_star=True),
                         None)

    def test_removable_gap_starts_here_does_not_detect_small_gaps(self):
        query_sequence = "fff........nnnnnn"
        hit_sequence   = "nnnnnnnnnnnnnnnnn"

        self.assertEqual(removable_gap_starts_here(query_sequence, hit_sequence, 3, gap_tolerated=7, always_remove_gap_with_star=True),
                         11)
        self.assertEqual(removable_gap_starts_here(query_sequence, hit_sequence, 3, gap_tolerated=10, always_remove_gap_with_star=True),
                         None)

    def test_removable_gap_starts_here_detects_any_gap_with_star_if_requested(self):
        query_sequence = "fff........nnnnnn"
        hit_sequence   = "nnnnnnn*nnnnnnnnn"

        self.assertEqual(removable_gap_starts_here(query_sequence, hit_sequence, 3, gap_tolerated=100, always_remove_gap_with_star=True),
                         11)
        self.assertEqual(removable_gap_starts_here(query_sequence, hit_sequence, 3, gap_tolerated=100, always_remove_gap_with_star=False),
                         None)

    def test_removable_gap_starts_here_works_with_gap_at_start(self):
        query_sequence = "........nnnnnn"
        hit_sequence   = "nnnnnnnnnnnnnn"
        self.assertEqual(removable_gap_starts_here(query_sequence, hit_sequence, 0, gap_tolerated=7, always_remove_gap_with_star=True),
                         8)

    def test_removable_gap_starts_here_works_with_gap_at_end(self):
        query_sequence = "nnnnnn........"
        hit_sequence   = "nnnnnnnnnnnnnn"
        self.assertEqual(removable_gap_starts_here(query_sequence, hit_sequence, 5, gap_tolerated=7, always_remove_gap_with_star=True),
                         None)
        self.assertEqual(removable_gap_starts_here(query_sequence, hit_sequence, 6, gap_tolerated=7, always_remove_gap_with_star=True),
                         len(query_sequence))

    def test_aligned_split_by_gap(self):
        x = Aligned(
            query_result_id = "q",
            query_length = 100,
            hit_id = "t",
            evalue = 0.001,
            evalue_cond = 0.0001,
            bitscore = 20,
            query_sequence = "fff.........nnn.........ii",
            query_positions_1b = [2,3,4,4,4,4,4,4,4,4,4,4,5,6,7,7,7,7,7,7,7,7,7,7,8,9],
            hit_sequence   = "nnnkkkkkkkkklllsssssssssii",
            hit_positions_1b = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26]
        )
        splitted = x.query_gap_removed()

        self.assertEqual(len(splitted), 3)
        self.assertEqual(splitted[0].query_sequence, "fff")
        self.assertEqual(splitted[0].query_positions_1b, [2, 3, 4])
        self.assertEqual(splitted[0].hit_sequence, "nnn")
        self.assertEqual(splitted[0].hit_positions_1b, [1, 2, 3])
        self.assertEqual(splitted[1].query_sequence, "nnn")
        self.assertEqual(splitted[1].query_positions_1b, [5, 6, 7])
        self.assertEqual(splitted[1].hit_sequence, "lll")
        self.assertEqual(splitted[1].hit_positions_1b, [13, 14, 15])
        self.assertEqual(splitted[2].query_sequence, "ii")
        self.assertEqual(splitted[2].query_positions_1b, [8, 9])
        self.assertEqual(splitted[2].hit_sequence, "ii")
        self.assertEqual(splitted[2].hit_positions_1b, [25, 26])
