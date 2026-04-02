from typing import Dict, List, Optional, Tuple
from Bio.Seq import Seq


def extract_subsequence(full_sequence: Optional[str], start_1_based: int, end_1_based: int) -> Optional[str]:
    if full_sequence is None:
        return None
    if start_1_based <= 0 or end_1_based <= 0:
        return None
    # coordinates are 1-based inclusive; order can be reversed depending on alignment direction
    left = min(start_1_based, end_1_based)
    right = max(start_1_based, end_1_based)
    if left > len(full_sequence):
        return None
    # slice is exclusive of end; adjust for 1-based inclusive
    return full_sequence[left - 1 : min(right, len(full_sequence))]


def extract_subsequence_strand_sensitive(full_sequence: Optional[str], start_1_based: int, end_1_based: int) -> Optional[str]:
    subs = extract_subsequence(full_sequence, start_1_based, end_1_based)
    if start_1_based > end_1_based:
        return str(Seq(subs).reverse_complement())
    return subs


def compute_three_frame_translations(full_seq, start, end):
    target_sequence = extract_subsequence_strand_sensitive(full_seq, start, end)
    if target_sequence is None:
        print("Cannot extract sequence using", len(full_seq), start, end)

    translations = []
    for frame in range(3):
        trim_right = (len(target_sequence)-frame)%3
        if trim_right > 0:
          frame_sequence = target_sequence[frame:-trim_right]
        else:
          frame_sequence = target_sequence[frame:]
        assert len(frame_sequence) % 3 == 0
        aa = Seq(frame_sequence).translate(to_stop=False) # Translate entire sequence, including stops

        if end > start:  # fwd strand
            translations.append((start+frame, end-trim_right, str(aa)))
        else:  # rev strand
            translations.append((start-frame, end+trim_right, str(aa)))

    return translations


def to_dna_coordinate(frame_dna_start, frame_dna_end, aa_from, aa_to):
    if frame_dna_end > frame_dna_start: # fwd strand
        dna_aa_from = frame_dna_start+(aa_from-1)*3
        dna_aa_to = frame_dna_start+aa_to*3-1
    else: # rev strand 
        dna_aa_from = frame_dna_start-(aa_from-1)*3
        dna_aa_to = frame_dna_start-aa_to*3+1
    return dna_aa_from, dna_aa_to


# pool multiple genome fasta files together, adding the genome accession to
# each contig name

import heapq

def partition_sequences(data_dict, target_l):
    # (length, key) sorted descending
    items = sorted(
        [(len(seq), key) for key, seq in data_dict.items()], 
        key=lambda x: x[0], 
        reverse=True
    )

    total_len = sum(item[0] for item in items)
    num_buckets = max(1, round(total_len / target_l))

    # (current_total_length, bucket_index, list_of_keys)
    buckets = [[0, i, []] for i in range(num_buckets)]
    heapq.heapify(buckets)

    for length, key in items:
        current_sum, idx, keys = heapq.heappop(buckets)
        keys.append(key)
        current_sum += length
        heapq.heappush(buckets, [current_sum, idx, keys])

    for i,b in enumerate(buckets):
        print(i, len(b[2]), b[0])
    return [b[2] for b in buckets]


def split_sequence_dictionary(data_dict, target_l):

    buckets = partition_sequences(data_dict, target_l)
    splitted = []

    for i, bucket in enumerate(buckets):
        new_dict = {x:data_dict[x] for x in bucket}
        splitted.append(new_dict)

    return splitted        



