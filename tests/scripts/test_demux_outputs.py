import os
import tempfile
import unittest
import subprocess
from pathlib import Path

from tangle.sequence import write_fasta_from_dict, read_fasta_as_dict
from tangle.detected import DetectedTable
from tangle.models import CSVSource


class TestDemuxScript(unittest.TestCase):

    def test_demux_script_creates_appropriate_protein_fastas(self):

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir = Path(temp_dir)
            in_tsv = temp_dir / "in.tsv"
            in_faa = temp_dir / "in.faa"
            out_tsv = temp_dir / "out.tsv"
            out_dir = temp_dir / "out"

            rows = [
                dict(detection_type="sequence", detection_method="hmm", batch="b1",
                     query_database="_", query_accession="g1|q1", query_type="contig",
                     target_database="_", target_accession="g1|t1", target_type="protein",
                     query_start=1, query_end=2, target_start=1, target_end=2),
                dict(detection_type="sequence", detection_method="hmm", batch="b1",
                     query_database="_", query_accession="g1|q2", query_type="contig",
                     target_database="_", target_accession="g1|t2", target_type="protein",
                     query_start=1, query_end=2, target_start=1, target_end=2),
                dict(detection_type="sequence", detection_method="hmm", batch="b1",
                     query_database="_", query_accession="g2|q3", query_type="contig",
                     target_database="_", target_accession="g2|t3", target_type="protein",
                     query_start=1, query_end=2, target_start=1, target_end=2),
                dict(detection_type="sequence", detection_method="hmm", batch="b1",
                     query_database="_", query_accession="g3|q5", query_type="contig",
                     target_database="_", target_accession="g3|t5", target_type="protein",
                     query_start=1, query_end=2, target_start=1, target_end=2),
                dict(detection_type="sequence", detection_method="hmm", batch="b1",
                     query_database="_", query_accession="g3|q6", query_type="contig",
                     target_database="_", target_accession="g3|t6", target_type="protein",
                     query_start=1, query_end=2, target_start=1, target_end=2)
            ]
            DetectedTable.write_tsv(str(in_tsv), rows)

            sequences = {
                "g1|t1": "A",
                "g1|t2": "S",
                "g2|t3": "L",
                "g3|t5": "M",
                "g3|t6": "Q"
            }
            write_fasta_from_dict(sequences, str(in_faa))

            cmd = ["python3", "scripts/demux-outputs.py", str(in_tsv), str(in_faa), str(out_tsv), str(out_dir)]

            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            self.assertEqual(result.stdout, "")

            g1_seq = read_fasta_as_dict((out_dir / "g1") / "proteins.faa")
            g2_seq = read_fasta_as_dict((out_dir / "g2") / "proteins.faa")
            g3_seq = read_fasta_as_dict((out_dir / "g3") / "proteins.faa")

            self.assertEqual(g1_seq, dict(t1="A", t2="S"))
            self.assertEqual(g2_seq, dict(t3="L"))
            self.assertEqual(g3_seq, dict(t5="M", t6="Q"))

            values = CSVSource(DetectedTable, str(out_tsv)).values()
            self.assertCountEqual([(row["target_database"], row["target_accession"]) for row in values],
                                  [("g1", "t1"), ("g1", "t2"), ("g2", "t3"), ("g3", "t5"), ("g3", "t6")])
