import os
import unittest

from needle.sequence import to_dna_coordinate

class TestProteinToDNACoordinateConverstion(unittest.TestCase):

    def test_to_dna_coordinate_converts_on_fwd_strand(self):
        # argments are frame_dna_start, frame_dna_end, aa_from, aa_to
        self.assertEqual(to_dna_coordinate(1, 100, 3, 6), (7, 18))
        self.assertEqual(to_dna_coordinate(2, 100, 3, 6), (8, 19))
        self.assertEqual(to_dna_coordinate(3, 100, 3, 6), (9, 20))
        self.assertEqual(to_dna_coordinate(4, 100, 3, 6), (10, 21))

    def test_to_dna_coordinate_converts_on_rev_strand(self):
        # argments are frame_dna_start, frame_dna_end, aa_from, aa_to
        self.assertEqual(to_dna_coordinate(100, 1, 3, 6), (94, 83))
        self.assertEqual(to_dna_coordinate( 99, 1, 3, 6), (93, 82))
        self.assertEqual(to_dna_coordinate( 98, 1, 3, 6), (92, 81))
        self.assertEqual(to_dna_coordinate( 97, 1, 3, 6), (91, 80))
