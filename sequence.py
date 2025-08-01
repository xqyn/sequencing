"""
project: sequencing
XQ - Leiden UMC
modifying func for dna sequence / dna list
update: 2025
    - march 14: add dna
    - april 22: compute_base_count
    - mei 2: hamming_dist
"""

from typing import NamedTuple, List, Dict
from itertools import zip_longest

# --- dna: modifying func for dna sequence
def ascii_to_phred(quality_string, offset=33):
    """Convert an ASCII quality to Phred scores"""
    return [ord(char) - offset for char in quality_string]


def complement(seq, reverse=False):
    """Return the complement of a sequence"""
    mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
               'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    complemented = ''.join([mapping.get(base, base) for base in seq])
    return complemented[::-1] if reverse else complemented


def rever_complement(seq):
    """Return the reverse complement of a sequence"""
    return complement(seq, reverse=True)


# --- hamming_dist
def hamming(ref_seq: str,
            comp_seq: str,
            position: bool = False,
            return_dict: bool = False) -> int:
    """
    Calculate the Hamming distance between two sequences.

    Parameters
    ----------
    ref_seq : str
        reference DNA sequence.
    comp_seq : str
        comparison DNA sequence.
    position : bool
        return the positions of the mismatches 
        (default: False).
    return_dict : bool
        returns a dict with key as Hamming distance, value as list of mismatch position
        (default: False).

    Returns
    -------
    int or list of int or dict
        number of positions (Hamming distance) between the two sequences.
        returns a list of positions where the mismatches occur.
        returns a dict with key as Hamming distance, value as list of mismatch position.

    Raises
    ------
    ValueError
        If sequences have different lengths or contain invalid nucleotides.
        If both position and return_dict are True.
        
    Examples
    --------
    hamming("ATCGATCGATCG", "ATGGATCGGTCG", position=True)
    [3, 9]
    hamming("ATCGATCGATCG", "ATGGATCGGTCG", return_dict=True)
    {2: [3, 9]}
    """
    
    if len(ref_seq) != len(comp_seq):
        raise ValueError("Sequences must have equal length for Hamming distance")
    
    valid_nucleotides = set('ACGTNQ-')
    if not (set(ref_seq).issubset(valid_nucleotides) and set(comp_seq).issubset(valid_nucleotides)):
        raise ValueError("Sequences must contain only A, C, G, T, or N")
    
    mismatches = [base for base, (ref, com) in enumerate(zip_longest(ref_seq, comp_seq), start=1) if ref != com]
    
    if return_dict:
        return {len(mismatches): mismatches}
    elif position:
        return mismatches
    
    return len(mismatches)


# --- compute_base_percent
def compute_base_count(seq_list: list,
                       bases: str = 'ATCGN-',
                       laplace: int = 0,
                       percentage: bool = False,
                       fraction: bool = False) -> dict:
    
    """
    calculate the count/percentage of each base at each base position 
    in a list of aligned DNA sequences.

    analyzes a list of DNA sequences and computes the number or proportion of 
    specified nucleotide bases (e.g., A, T, C, G, N, -) at each position. 
    All sequences are assumed to be aligned and of the same length.
    
    Parameters
    ----------
    seq_list : list
        A list of DNA sequences string, equal length.
    bases : str
        string of valid dna base.
        (default: "ATCGN-").
    laplace : str
        Laplace’s Rule of Succession: add a pseudocount into the matrix
        (default: 0).
    percentage : bool
        whether to returns the count or percentage occurrence
        (default: False)
    fraction : bool
        Return fractions (sum to 1) instead of percentages 
        (default: False).

    Returns
    -------
    dict of list of int or float
        a dictionary where each key is a base character from `bases`, 
        and the corresponding value is a list of integers (counts) or floats (percentages)
        representing the base frequency at each position across all sequences.
    
    Examples
    --------
    seq_list = ["ATCGATCGATCG",
                "ATGGAATG-GC",
                "ATCGATCGGTC",
                "ATG-TGAT-GC",
                "ATCGNNAT-NG"]
    compute_base_count(seq_list)
    compute_base_count(seq_list, percentage=True)
    compute_base_count(seq_list, laplace = 1, percentage=True)
    """
    if not seq_list:
        return {base: [] for base in bases}
    
    # assume length of sequence in seq_list are equal
    length = len(seq_list[0])
    # init dict for base count for each position
    # initialize with Laplace pseudocount
    base_dict = {base: [laplace] * length for base in bases}
    
    # count occurences of each base from bases at each position
    for seq in seq_list:
        for pos, base in enumerate(seq):
            if base in bases:
                base_dict[base][pos] += 1
    
    if percentage or fraction:
        total_counts = len(seq_list) + (laplace * len(bases))
        multiplier = 1 if fraction else 100
        return {
            base: [round(count / total_counts * multiplier, 5) for count in counts]
            for base, counts in base_dict.items()
        }
    return base_dict
