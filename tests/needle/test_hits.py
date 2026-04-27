import os
import io
import tempfile
import unittest

from needle.match import ProteinHit, Match
from needle.hits import write_fasta_record, export_protein_hits


class TestExport(unittest.TestCase):

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
            export_protein_hits("GENOME1", [pm], p1, p2, append=True)
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
