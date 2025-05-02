#!/bin/bash python
"""
project: sequencing
XQ - Leiden UMC
motif algorithm
update: 2025
    - march 19: pattern_in_seq
    - mei 2: motif_score

Terminologies:
    - motif/pattern:    a string of base, usually short (e.g: k-mer)
    - dna/sequence:     a string of base, longer which may contrain motifs inside (e.g: dna sequence)
    - motifs:           a collection of motifs (e.g: list of motif [motif_1, motif_2,... motif_n])
"""

from collections import Counter
from seq import hamming


# --- pattern_in_seq
def pattern_in_seq(pattern: str,
                   sequence: str,
                   upper: bool = True) -> tuple[int, Counter]:
    """
    find motifs in the sequence with the lowest Hamming distance to the given pattern.

    Parameters
    ----------
    pattern : str
        pattern to match against.
    sequence : str
        string containing potential motifs.
    upper : bool
        If True, converts pattern and sequence to uppercase before processing. Default is True.

    Returns
    -------
    dict
        Dictionary with minimum Hamming distance as key and a dictionary of motifs with that
        distance as value.

    Raises
    ------
    IndexError
        If pattern length is greater than sequence length.

    Examples
    --------
    >>> pattern_in_seq("ATG", "ATGATGCGAT")
    {0: {'ATG': 2}}
    >>> pattern_in_seq("ATG", "ATCGTGCGAT")
    {1: {'ATC': 1, 'GTG': 1}}
    """
    if upper:
        pattern = pattern.upper()
        sequence = sequence.upper()
      
    ptn_len = len(pattern)
    min_dist = ptn_len
    k_mer_counts = {}
    
    # Generate k-mers and track those with minimum distance
    for idx in range(len(sequence) - ptn_len + 1):
        k_mer = sequence[idx:idx + ptn_len]
        dist = hamming(pattern, k_mer)
        if dist <= min_dist:
            if dist < min_dist:
                k_mer_counts = {}
                min_dist = dist
            k_mer_counts[k_mer] = k_mer_counts.get(k_mer, 0) + 1
    
    return {min_dist: k_mer_counts}


# --- motif_score
def motif_score(motifs: list,
                upper: bool = True,
                col: bool = True) -> int:
    """
    calculate the sum of non-consensus bases for motifs list.

    Parameters
    ----------
    motifs : list
        list of motif strings (equal-length DNA sequences).
    upper : bool
        converts all bases to uppercase before scoring. 
        (default: True).
    col : bool
       score by column (most common base per column); if False, score by row
        (non-matching bases per sequence). 
        (default is True).

    Returns
    -------
    int
        total score, summing the number of bases that DIFFER from the most common
        base in each position.

    Raises
    ------
    IndexError
        If motifs is empty or contains sequences of unequal length.

    Examples
    --------
    >>> motifs = [
                "TCGGGGgTTTtt",
                "cCGGtGAcTTaC",
                "aCGGGGATTTtC",
                "TtGGGGAcTTtt",
                "aaGGGGAcTTCC",
                "TtGGGGAcTTCC",
                "TCGGGGATTcat",
                "TCGGGGATTcCt",
                "TaGGGGAacTaC",
                "TCGGGtATaaCC"
                ]
    >>> motif_score(motifs)
    >>> motif_score(motifs, col=False)
    """
    if upper:
        motifs = [motif.upper() for motif in motifs]
    
    k = len(motifs[0])
    score = 0
    if col:
        for i in range(k):
            column = [motif[i] for motif in motifs]
            most_common = max(set(column), key=column.count)
            score += sum(1 for base in column if base != most_common)
    else:
        for motif in motifs:
            for i in range(k):
                column = [m[i] for m in motifs]
                most_common = max(set(column), key=column.count)
                if motif[i] != most_common:
                    score += 1
    return score