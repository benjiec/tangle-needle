import os
import tempfile
import unittest
from needle.detect import Results, extract_fragments, get_aa_sequences
import needle.detect as detect_mod


class TestParseDetectResults(unittest.TestCase):

    def test_parse_ncbi_header_and_extract_sequences_reverse(self):
        """Verify reverse-direction hit: target sequence normalized to 5'->3' and reverse-complemented for translation."""

        # Create synthetic FASTA and TSV with NCBI-style headers; reverse-direction on target (sstart > send)
        with tempfile.TemporaryDirectory() as tmpdir:
            target_fasta_path = os.path.join(tmpdir, "target.faa")
            results_tsv_path = os.path.join(tmpdir, "results.tsv")

            # >S1 DNA: AAACCCGGGTTT
            with open(target_fasta_path, "w") as f:
                f.write(">S1 synthetic subject\n")
                f.write("AAACCCGGGTTT\n")

            header = "\t".join(Results.PRODUCER_HEADER)
            # qseqid sseqid evalue pident qstart qend sstart send sseq
            row = ["Q1", "S1", "1e-5", "99.9", "2", "5", "10", "7", "ABCD"]
            with open(results_tsv_path, "w") as f:
                f.write(header + "\n")
                f.write("\t".join(row) + "\n")

            res = Results(results_tsv_path, target_fasta_path)
            matches = res.matches()
            self.assertEqual(len(matches), 1)
            m = matches[0]
            self.assertEqual(m.query_accession, "Q1")
            self.assertEqual(m.target_accession, "S1")
            self.assertAlmostEqual(m.e_value, 1e-5, places=12)
            self.assertAlmostEqual(m.identity, 99.9, places=3)
            self.assertEqual(m.query_start, 2)
            self.assertEqual(m.query_end, 5)
            self.assertEqual(m.target_start, 10)
            self.assertEqual(m.target_end, 7)

            # Target 7..10 on AAACCCGGGTTT -> "GGGT"; reverse-complement -> "ACCC"
            self.assertEqual(m.target_sequence, "ACCC")
            # Ensure translated target (revcomp path) is a valid amino-acid sequence length <= query segment
            _ = m.target_sequence_translated()
            # Matched sequence preserved
            self.assertEqual(m.matched_sequence, "ABCD")

    def test_parse_ncbi_header_and_extract_sequences_forward(self):
        """Verify forward-direction hit: target sequence 5'->3' matches the query AA after translation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            target_fasta_path = os.path.join(tmpdir, "t.fna")
            results_tsv_path = os.path.join(tmpdir, "r.tsv")

            # Target DNA forward: ATG GAA TTT -> MEF at positions 5..13
            with open(target_fasta_path, "w") as f:
                f.write(">Tfwd\n")
                f.write("NNNNATGGAATTTNNNN\n")  # 1..4 N; 5..13 coding; 14..17 N

            header = "\t".join(Results.PRODUCER_HEADER)
            # Forward direction: sstart < send
            row = ["Qfwd", "Tfwd", "0", "100.0", "1", "3", "5", "13", "MEF"]
            with open(results_tsv_path, "w") as f:
                f.write(header + "\n")
                f.write("\t".join(row) + "\n")

            res = Results(results_tsv_path, target_fasta_path)
            ms = res.matches()
            self.assertEqual(len(ms), 1)
            m = ms[0]
            # Target DNA 5..13 5'->3'
            self.assertEqual(m.target_sequence, "ATGGAATTT")
            # Translated target equals matched_sequence equals query segment
            self.assertEqual(m.target_sequence_translated(), "MEF")
            self.assertEqual(m.matched_sequence, "MEF")


class TestHMMSearchGenome(unittest.TestCase):

    def setUp(self):
        self.orig = detect_mod.hmmsearch_sequence_dict
        self.hmmsearch_calls_translated = []

        def _fake_hmmsearch(hmm_file, translated): 
            self.hmmsearch_calls_translated.append(translated)
            return [
                dict(target_name = k,
                     target_accession = "t1",
                     query_name = "q",
                     query_accession = "q1",
                     seq_evalue = 0.001,
                     seq_score = 20,
                     dom_evalue = 0.001,
                     dom_score = 20,
                     query_length = 60,
                     hmm_from = 3,
                     hmm_to = 6,
                     target_length = 7,
                     ali_from = 2,
                     ali_to = 5 
                ) for k,v in translated.items()
            ]

        detect_mod.hmmsearch_sequence_dict = _fake_hmmsearch

    def tearDown(self):
        detect_mod.hmmsearch_sequence_dict = self.orig

    def test_extract_fragments_splits_by_star_and_reports_dna_coordinates_correctly_fwd_strand(self):

        seq = "F"*10+"*"+"G"*5+"*"+"L"
        fragments = extract_fragments("t1", 123, 345, seq)
        self.assertEqual(len(fragments), 3)
        self.assertEqual(fragments[0][0], "t1")
        self.assertEqual(fragments[0][1], 123)
        self.assertEqual(fragments[0][2], 123+10*3-1)
        self.assertEqual(fragments[0][3], "F"*10)
        self.assertEqual(fragments[1][0], "t1")
        self.assertEqual(fragments[1][1], 123+11*3)
        self.assertEqual(fragments[1][2], 123+(11+5)*3-1)
        self.assertEqual(fragments[1][3], "G"*5)
        self.assertEqual(fragments[2][0], "t1")
        self.assertEqual(fragments[2][1], 123+(11+6)*3)
        self.assertEqual(fragments[2][2], 123+(11+6+1)*3-1)
        self.assertEqual(fragments[2][3], "L")

    def test_extract_fragments_splits_handles_starting_star(self):

        seq = "*"+"F"*10+"*"+"G"*5+"*"+"L"
        fragments = extract_fragments("t1", 120, 345, seq)
        self.assertEqual(len(fragments), 3)
        self.assertEqual(fragments[0][0], "t1")
        self.assertEqual(fragments[0][1], 123)
        self.assertEqual(fragments[0][2], 123+10*3-1)
        self.assertEqual(fragments[0][3], "F"*10)
        self.assertEqual(fragments[1][0], "t1")
        self.assertEqual(fragments[1][1], 123+11*3)
        self.assertEqual(fragments[1][2], 123+(11+5)*3-1)
        self.assertEqual(fragments[1][3], "G"*5)
        self.assertEqual(fragments[2][0], "t1")
        self.assertEqual(fragments[2][1], 123+(11+6)*3)
        self.assertEqual(fragments[2][2], 123+(11+6+1)*3-1)
        self.assertEqual(fragments[2][3], "L")

    def test_get_aa_sequences_returns_fragments_using_six_frames(self):

        seq = "ATAA"*4  # this will have a lot of TAAs in some frames
        # fwd: ATAAATAAATAAATAA  1..16
        #  f0: ataAATaaaTAAata   => 2
        #  f1:  taaATAaatAAAtaa  => 1
        #  f2:   aaaTAAataAAT    => 2
        # rev: TTATTTATTTATTTAT  16..1
        #  r0: ttaTTTattTATtta   => 1
        #  r1:  TATttaTTTattTAT  => 1
        #  r2:   ATTtatTTAttt    => 1

        fragments = get_aa_sequences("t1", seq, min_aa_length=1)

        self.assertEqual(len(fragments), 2+1+2+1+1+1)
        # f0.1
        self.assertEqual(fragments[0][1], 1)
        self.assertEqual(fragments[0][2], 9)
        # f0.2
        self.assertEqual(fragments[1][1], 13)
        self.assertEqual(fragments[1][2], 15)
        # f1.1
        self.assertEqual(fragments[2][1], 5)
        self.assertEqual(fragments[2][2], 13)
        # f2.1
        self.assertEqual(fragments[3][1], 3)
        self.assertEqual(fragments[3][2], 5)
        # f2.2
        self.assertEqual(fragments[4][1], 9)
        self.assertEqual(fragments[4][2], 14)
        # r0.1
        self.assertEqual(fragments[5][1], 16)
        self.assertEqual(fragments[5][2], 2)
        # r1.1
        self.assertEqual(fragments[6][1], 15)
        self.assertEqual(fragments[6][2], 1)
        # r2.1
        self.assertEqual(fragments[7][1], 14)
        self.assertEqual(fragments[7][2], 3)

    def test_get_aa_sequences_returns_fragments_correctly_even_when_windowing(self):

        seq = "ATAA"*4  # this will have a lot of TAAs in some frames
        fragments = get_aa_sequences("t1", seq, win=8, win_overlap=4, min_aa_length=1)

        # first window, 1..12
        # fwd: ATAAATAAATAA 1..12
        #  f0: ataAATaaaTAA      => 1
        #  f1:  taaATAaat        => 1
        #  f2:   aaaTAAata       => 2
        # rev: TTATTTATTTAT 12..1
        #  r0: ttaTTTattTAT      => 1
        #  r1:  TATttaTTT        => 1
        #  r2:   ATTtatTTA       => 1

        # second window, 9..16
        # fwd: ATAAATAA     9..16
        #  f0: ataAAT            => 1
        #  f1:  taaATA           => 1
        #  f2:   aaaTAA          => 1
        # rev: TTATTTAT     16..9
        #  r0: ttaTTT            => 1
        #  r1:  TATtta           => 1
        #  r2:   ATTtat          => 1

        self.assertEqual(len(fragments), 1+1+2+ 1+1+1+ 1+1+1+ 1+1+1)
        # w1 f0.1
        self.assertEqual(fragments[0][1], 1)
        self.assertEqual(fragments[0][2], 9)
        # w1 f1.1
        self.assertEqual(fragments[1][1], 5)
        self.assertEqual(fragments[1][2], 10)
        # w1 r0.1
        self.assertEqual(fragments[4][1], 12)
        self.assertEqual(fragments[4][2], 1)
        # w2 f0.1
        self.assertEqual(fragments[7][1], 9)
        self.assertEqual(fragments[7][2], 14)
        # w2 r0.1
        self.assertEqual(fragments[10][1], 16)
        self.assertEqual(fragments[10][2], 11)
