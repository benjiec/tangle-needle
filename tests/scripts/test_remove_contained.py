import os
import tempfile
import unittest
import subprocess
from pathlib import Path

from tangle.sequence import write_fasta_from_dict, read_fasta_as_dict
from tangle.detected import DetectedTable
from tangle.models import CSVSource


class TestRemoveContained(unittest.TestCase):

    def test_remove_contained_removes_contained_sequence_at_same_locus(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir = Path(temp_dir)
            in_tsv = temp_dir / "proteins.tsv"
            in_faa = temp_dir / "proteins.faa"

            rows = [
                dict(detection_type="sequence", detection_method="hmm", batch="b1",
                     query_database="g1", query_accession="q1", query_type="contig",
                     target_database="g1", target_accession="t1", target_type="protein",
                     query_start=10, query_end=19, target_start=1, target_end=2),
                dict(detection_type="sequence", detection_method="hmm", batch="b1",
                     query_database="g1", query_accession="q1", query_type="contig",
                     target_database="g1", target_accession="t2", target_type="protein",
                     query_start=5, query_end=23, target_start=1, target_end=2),
                dict(detection_type="sequence", detection_method="hmm", batch="b1",
                     query_database="g3", query_accession="q5", query_type="contig",
                     target_database="g3", target_accession="t5", target_type="protein",
                     query_start=10, query_end=19, target_start=1, target_end=2),
                dict(detection_type="sequence", detection_method="hmm", batch="b1",
                     query_database="g3", query_accession="q5", query_type="contig",
                     target_database="g3", target_accession="t6", target_type="protein",
                     query_start=1, query_end=15, target_start=1, target_end=2)
            ]

            # first two, same locus
            assert rows[0]["query_accession"] == rows[1]["query_accession"]
            assert rows[0]["query_start"] > rows[1]["query_start"]
            assert rows[0]["query_end"] < rows[1]["query_end"]
            # second two, overlapping loci (i.e. not contained)
            assert rows[2]["query_accession"] == rows[3]["query_accession"]
            assert rows[2]["query_start"] > rows[3]["query_start"]
            assert rows[2]["query_end"] > rows[3]["query_end"]  # not contained!

            DetectedTable.write_tsv(str(in_tsv), rows)

            sequences = {
                "t1": "AAA",
                "t2": "SAAAS",
                "t5": "MMM",
                "t6": "QMMMQ"
            }
            write_fasta_from_dict(sequences, str(in_faa))

            cmd = ["python3", "scripts/remove-contained.py", str(temp_dir)]
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)

            values = CSVSource(DetectedTable, str(in_tsv)).values()
            self.assertCountEqual([(row["target_database"], row["target_accession"]) for row in values],
                                  [("g1", "t2"), ("g3", "t5"), ("g3", "t6")])

            updated_seq = read_fasta_as_dict(str(in_faa))
            self.assertEqual(updated_seq, dict(t2="SAAAS", t5="MMM", t6="QMMMQ"))

    def test_remove_contained_works_with_reversed_query_coordinates(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir = Path(temp_dir)
            in_tsv = temp_dir / "proteins.tsv"
            in_faa = temp_dir / "proteins.faa"

            rows = [
                dict(detection_type="sequence", detection_method="hmm", batch="b1",
                     query_database="g1", query_accession="q1", query_type="contig",
                     target_database="g1", target_accession="t1", target_type="protein",
                     query_start=19, query_end=10, target_start=1, target_end=2),
                dict(detection_type="sequence", detection_method="hmm", batch="b1",
                     query_database="g1", query_accession="q1", query_type="contig",
                     target_database="g1", target_accession="t2", target_type="protein",
                     query_start=23, query_end=5, target_start=1, target_end=2),
            ]

            DetectedTable.write_tsv(str(in_tsv), rows)
            sequences = { "t1": "AAA", "t2": "SAAAS" }
            write_fasta_from_dict(sequences, str(in_faa))

            cmd = ["python3", "scripts/remove-contained.py", str(temp_dir)]
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)

            values = CSVSource(DetectedTable, str(in_tsv)).values()
            self.assertCountEqual([(row["target_database"], row["target_accession"]) for row in values], [("g1", "t2")])

            updated_seq = read_fasta_as_dict(str(in_faa))
            self.assertEqual(updated_seq, dict(t2="SAAAS"))

            # make only one of them reversed
            t = rows[0]["query_end"]
            rows[0]["query_end"] = rows[0]["query_start"]
            rows[0]["query_start"] = t

            DetectedTable.write_tsv(str(in_tsv), rows)
            sequences = { "t1": "AAA", "t2": "SAAAS" }
            write_fasta_from_dict(sequences, str(in_faa))

            cmd = ["python3", "scripts/remove-contained.py", str(temp_dir)]
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)

            values = CSVSource(DetectedTable, str(in_tsv)).values()
            self.assertCountEqual([(row["target_database"], row["target_accession"]) for row in values], [("g1", "t1"), ("g1", "t2")])

            updated_seq = read_fasta_as_dict(str(in_faa))
            self.assertEqual(updated_seq, dict(t1="AAA", t2="SAAAS"))

    def test_remove_contained_keeps_contained_from_very_diff_locus(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir = Path(temp_dir)
            in_tsv = temp_dir / "proteins.tsv"
            in_faa = temp_dir / "proteins.faa"

            rows = [
                dict(detection_type="sequence", detection_method="hmm", batch="b1",
                     query_database="g1", query_accession="q1", query_type="contig",
                     target_database="g1", target_accession="t1", target_type="protein",
                     query_start=19, query_end=10, target_start=1, target_end=2),
                dict(detection_type="sequence", detection_method="hmm", batch="b1",
                     query_database="g1", query_accession="q1", query_type="contig",
                     target_database="g1", target_accession="t2", target_type="protein",
                     query_start=23, query_end=5, target_start=1, target_end=2),
            ]

            DetectedTable.write_tsv(str(in_tsv), rows)
            sequences = { "t1": "AAA", "t2": "SAAAS" }
            write_fasta_from_dict(sequences, str(in_faa))

            cmd = ["python3", "scripts/remove-contained.py", str(temp_dir)]
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)

            values = CSVSource(DetectedTable, str(in_tsv)).values()
            self.assertCountEqual([(row["target_database"], row["target_accession"]) for row in values], [("g1", "t2")])

            updated_seq = read_fasta_as_dict(str(in_faa))
            self.assertEqual(updated_seq, dict(t2="SAAAS"))

            # move location dramatically
            rows[0]["query_start"] += 100
            rows[0]["query_end"] += 100

            DetectedTable.write_tsv(str(in_tsv), rows)
            sequences = { "t1": "AAA", "t2": "SAAAS" }
            write_fasta_from_dict(sequences, str(in_faa))

            cmd = ["python3", "scripts/remove-contained.py", str(temp_dir)]
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)

            values = CSVSource(DetectedTable, str(in_tsv)).values()
            self.assertCountEqual([(row["target_database"], row["target_accession"]) for row in values], [("g1", "t1"), ("g1", "t2")])

            updated_seq = read_fasta_as_dict(str(in_faa))
            self.assertEqual(updated_seq, dict(t1="AAA", t2="SAAAS"))
