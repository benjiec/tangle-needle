from typing import Dict, List, Optional, Tuple
from tangle.sequence import extract_subsequence, extract_subsequence_strand_sensitive, compute_three_frame_translations


def uniq_aa(s):
    lookup = {c:1 for c in s}
    return len(lookup.keys())

def is_aa_sequence_low_complexity(s):
    return uniq_aa(s) <= 8


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



