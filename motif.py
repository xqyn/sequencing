#!/bin/bash python

'''
Motif finding - XQ Ng - 19 March 2025
Leiden UMC
''' 

from collections import Counter
from itertools import zip_longest

# --------------------------------------------------
def hamming_dist(ref_seq: str, 
                comp_seq: str) -> int:
    """ Calculate Hamming distance 
    
    params:
        ref_seq: reference sequence
        comp_seq: comparison sequence
    
    return
        hamming distace: distance between reference and comparison sequence
    """
    
    return sum(map(lambda base: base[0] != base[1], zip_longest(ref_seq, comp_seq)))


def motif_finding(pattern: str, 
                  sequence: str,
                  to_upper: bool = True) -> tuple[int, Counter]:
    """ Find motifs in the sequence with the lowest Hamming distance to the given pattern.
    Suppose that pattern <= sequence in length
    
    params:
        pattern: pattern to match against
        sequence: string containing potential motifs
    
    return
        (minimum Hamming distance, Counter of motifs with minimum distance)
    """
    if to_upper:
        pattern = pattern.upper()
        sequence = sequence.upper()
    
    ptn_len = len(pattern)
    min_dist = ptn_len
    k_mer_candidates = []
    
    # Generate k-mers and track those with minimum distance
    for idx in range(len(sequence) - ptn_len + 1):
        k_mer = sequence[idx:idx + ptn_len]
        dist = hamming_dist(pattern, k_mer)
        if dist < min_dist:
            min_dist = dist
            k_mer_candidates = [k_mer]
        elif dist == min_dist:
            k_mer_candidates.append(k_mer)
    
    return min_dist, Counter(k_mer_candidates)
