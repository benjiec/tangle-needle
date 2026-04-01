import re
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from typing import List
from Bio import SearchIO


def run_command(cmd: str):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def hmmfetch(hmm_file_name, acc, output_file_name):
    cmd = ["hmmfetch", "-o", output_file_name, hmm_file_name, acc]
    run_command(cmd)


class HMMCollection(object):

    def __init__(self, big_hmm_file, accession_ids):
        self.__by_accession = {}

        accession_ids = list(set(accession_ids))
        with tempfile.TemporaryDirectory(delete=False) as temp_dir:
            self.__temp_dir = temp_dir
            for acc in accession_ids:
                with tempfile.NamedTemporaryFile(dir=temp_dir, delete=False, suffix=".hmm") as tmpf:
                    tmpf.close()
                    # print("fetching hmm into temp file", tmpf.name)
                    try: # file may not be there
                        hmmfetch(big_hmm_file, acc, tmpf.name)
                        self.__by_accession[acc] = tmpf.name
                    except:
                        os.remove(tmpf.name)

    def get(self, accession):
        return self.__by_accession[accession] if accession in self.__by_accession else None

    def clean(self):
        for fn in self.__by_accession.values():
            # print("removing", fn)
            os.remove(fn)
        # print("removing", self.__temp_dir)
        shutil.rmtree(self.__temp_dir)
        self.__by_accession = None
        self.__temp_dir = None


def parse_hmmsearch_domtbl(domtbl_path):
    expected_header = "# target name  accession  tlen  query name  accession  qlen  E-value  score  bias  #  of  c-Evalue  i-Evalue  score  bias  from  to  from  to"
    idx_target = 0
    idx_target_acc = 1
    idx_t_len = 2
    idx_query = 3
    idx_query_acc = 4
    idx_q_len = 5
    idx_seq_eval = 6
    idx_seq_score = 7
    idx_dom_eval_cond = 11
    idx_dom_eval = 12
    idx_dom_score = 13
    idx_h_from = 15
    idx_h_to = 16
    idx_a_from = 17
    idx_a_to = 18

    # first, sanity check these indices
    expected_header_parts = re.split(r'\s\s+', expected_header)
    assert expected_header_parts[idx_target] == "# target name"
    assert expected_header_parts[idx_target_acc] == "accession"
    assert expected_header_parts[idx_t_len] == "tlen"
    assert expected_header_parts[idx_query] == "query name"
    assert expected_header_parts[idx_query_acc] == "accession"
    assert expected_header_parts[idx_q_len] == "qlen"
    assert expected_header_parts[idx_seq_eval] == "E-value"
    assert expected_header_parts[idx_seq_score] == "score"
    assert expected_header_parts[idx_dom_eval_cond] == "c-Evalue"
    assert expected_header_parts[idx_dom_eval] == "i-Evalue"
    assert expected_header_parts[idx_dom_score] == "score"
    assert expected_header_parts[idx_h_from] == "from"
    assert expected_header_parts[idx_h_to] == "to"
    assert expected_header_parts[idx_a_from] == "from"
    assert expected_header_parts[idx_a_to] == "to"

    has_headers = False
    expected_header = " ".join(expected_header.split())

    matches = []
    with open(domtbl_path, "r") as domf:
        for line in domf:
            if " ".join(line.split()).startswith(expected_header):
                has_headers = True
            if not line or line.startswith("#") or has_headers is False:
                continue
            parts = line.strip().split()
            match = dict(
                target_name = parts[idx_target],
                target_accession = parts[idx_target_acc],
                query_name = parts[idx_query],
                query_accession = parts[idx_query_acc],
                seq_evalue = float(parts[idx_seq_eval]),
                seq_score = float(parts[idx_seq_score]),
                dom_evalue = float(parts[idx_dom_eval]),
                dom_evalue_cond = float(parts[idx_dom_eval_cond]),
                dom_score = float(parts[idx_dom_score]),
                query_length = int(parts[idx_q_len]),
                hmm_from = int(parts[idx_h_from]),
                hmm_to = int(parts[idx_h_to]),
                target_length = int(parts[idx_t_len]),
                ali_from = int(parts[idx_a_from]),
                ali_to = int(parts[idx_a_to])
            )
            matches.append(match)

    assert has_headers
    return matches


def removable_gap_starts_here(query_sequence, hit_sequence, pos, gap_tolerated, always_remove_gap_with_star):
    gap_started = pos
    seen_star_on_hit = False
    for i in range(pos, len(query_sequence)):
        if query_sequence[i] in ".-":
            if always_remove_gap_with_star and hit_sequence[i] == "*":
                seen_star_on_hit = True
        else:
            if i-gap_started > gap_tolerated or \
               (i-gap_started > 1 and seen_star_on_hit):
                return i
            return None
    return len(query_sequence)


@dataclass
class Aligned(object):
    query_result_id: str
    query_length: int
    hit_id: str
    evalue: float
    evalue_cond: float
    bitscore: float
    query_sequence: List[str]
    query_positions_1b: List[int]
    hit_sequence: List[str]
    hit_positions_1b: List[int]

    def match_row(self):
        return dict(
            target_name = self.hit_id,
            target_accession = self.hit_id,
            query_name = self.query_result_id,
            query_accession = self.query_result_id,
            seq_evalue = None,
            seq_score = None,
            dom_evalue = self.evalue,
            dom_evalue_cond = self.evalue_cond,
            dom_score = self.bitscore,
            query_length = self.query_length,
            hmm_from = self.query_positions_1b[0],
            hmm_to = self.query_positions_1b[-1],
            target_length = None,
            ali_from = self.hit_positions_1b[0],
            ali_to = self.hit_positions_1b[-1],
        )

    @staticmethod
    def from_hsp(query_result_id, query_length, hit_id, hsp):
        assert len(hsp.query.seq) == len(hsp.hit.seq)

        query_sequence = []
        hit_sequence = []
        query_positions_1b = []
        hit_positions_1b = []

        if hsp.query_start <= hsp.query_end:
            _query_incr = 1
        else:
            _query_incr = -1
        if hsp.hit_start <= hsp.hit_end:
            _hit_incr = 1
        else:
            _hit_incr = -1

        # we are recording 1-based positions, BioPython convention is [start..end) and 0 based
        next_query_pos = hsp.query_start+1
        next_hit_pos = hsp.hit_start+1

        for i,(q,t) in enumerate(zip(str(hsp.query.seq), str(hsp.hit.seq))):
           query_sequence.append(q)
           hit_sequence.append(t)
           query_positions_1b.append(next_query_pos)
           hit_positions_1b.append(next_hit_pos)

           if q in "-.":  # there is an insertion on target sequence
               pass       # query position does not advance
           else:
               next_query_pos += _query_incr

           if t in "-.":  # there is a deletion on target sequence
               pass       # target position does not advance
           else:
               next_hit_pos += _hit_incr

        # comparing 1based with 0based+1
        assert query_positions_1b[-1] == hsp.query_end
        assert hit_positions_1b[-1] == hsp.hit_end

        return Aligned(
            query_result_id = query_result_id,
            query_length = query_length,
            hit_id = hit_id,
            evalue = hsp.evalue,
            evalue_cond = hsp.evalue_cond,
            bitscore = hsp.bitscore,
            query_sequence = query_sequence,
            query_positions_1b = query_positions_1b,
            hit_sequence = hit_sequence,
            hit_positions_1b = hit_positions_1b
        )

    def query_gap_removed(self, gap_tolerated=8, always_remove_gap_with_star=True):

        aligned = []

        last_boundary_i = 0
        i = 0
        while i < len(self.query_sequence):
            pos_after_gap = removable_gap_starts_here(self.query_sequence, self.hit_sequence, i, gap_tolerated, always_remove_gap_with_star)
            if pos_after_gap is not None:
                aligned.append(Aligned(
                    query_result_id = self.query_result_id,
                    query_length = self.query_length,
                    hit_id = self.hit_id,
                    evalue = self.evalue,
                    evalue_cond = self.evalue_cond,
                    bitscore = self.bitscore,
                    query_sequence = self.query_sequence[last_boundary_i:i],
                    query_positions_1b = self.query_positions_1b[last_boundary_i:i],
                    hit_sequence = self.hit_sequence[last_boundary_i:i],
                    hit_positions_1b = self.hit_positions_1b[last_boundary_i:i]
                ))
                last_boundary_i = pos_after_gap
                i = pos_after_gap
            else:
                i += 1

        if i > last_boundary_i:
            aligned.append(Aligned(
                query_result_id = self.query_result_id,
                query_length = self.query_length,
                hit_id = self.hit_id,
                evalue = self.evalue,
                evalue_cond = self.evalue_cond,
                bitscore = self.bitscore,
                query_sequence = self.query_sequence[last_boundary_i:i],
                query_positions_1b = self.query_positions_1b[last_boundary_i:i],
                hit_sequence = self.hit_sequence[last_boundary_i:i],
                hit_positions_1b = self.hit_positions_1b[last_boundary_i:i]
            ))

        return aligned


def parse_hmmsearch_output(output_path):
    alignments = []
    with open(output_path, "r") as handle:
        for query_result in SearchIO.parse(handle, "hmmer3-text"):
            for hit in query_result.hits:
                for hsp in hit.hsps:
                    full_alignment = Aligned.from_hsp(query_result.id, query_result.seq_len, hit.id, hsp)
                    alignments.extend(full_alignment.query_gap_removed())
    return [x.match_row() for x in alignments]


def hmmsearch_file(hmm_file_name, fasta_path, cutoff=False, gap_removal=True, cpu=None):
    with tempfile.NamedTemporaryFile(delete=False, suffix=".domtbl", mode="w") as domtbl_f:
        domtbl_f.close()
        with tempfile.NamedTemporaryFile(delete=False, suffix=".txt", mode="w") as out_f:
            out_f.close()
            cmd = ["hmmsearch"]
            if cpu is not None:
                cmd.extend(["--cpu", str(cpu)])
            if cutoff:
                cmd.append("--cut_ga")
            cmd.extend(["-o", out_f.name, "--domtblout", domtbl_f.name, hmm_file_name, fasta_path])

            run_command(cmd)

            if gap_removal:
		# parse the alignments and return alignments w/o gaps on query,
		# as those may be introns
                res = parse_hmmsearch_output(out_f.name)
            else:
                res = parse_hmmsearch_domtbl(domtbl_f.name)

            os.remove(domtbl_f.name)
            os.remove(out_f.name)
            return res


def hmmsearch(hmm_file_name, sequences, cutoff=False, gap_removal=True, cpu=None):
    if len(sequences) == 0:
        return []
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, "cands.faa")
        with open(fasta_path, "w") as f:
            for i, cand in enumerate(sequences):
                f.write(f">cand_{i}\n{cand}\n")
        return hmmsearch_file(hmm_file_name, fasta_path, cutoff=cutoff, gap_removal=gap_removal, cpu=cpu)


def hmmsearch_sequence_dict(hmm_file_name, fasta_dict, cutoff=False, gap_removal=True, cpu=None):
    if len(fasta_dict.keys()) == 0:
        return []
    with tempfile.NamedTemporaryFile(delete=False, suffix=".faa", mode="w") as tmpf:
        for acc, sequence in fasta_dict.items():
            tmpf.write(f">{acc}\n{sequence}\n")
        tmpf.close()
        res = hmmsearch_file(hmm_file_name, tmpf.name, cutoff=cutoff, gap_removal=gap_removal, cpu=cpu)
        os.remove(tmpf.name)
        return res
