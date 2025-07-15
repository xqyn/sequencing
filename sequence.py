"""
project: sequencing
XQ - Leiden UMC
modifying func for dna sequence / dna list
update: 2025
    - april 22: compute_base_count
    - mei 2: hamming_dist
    - juli 15: cigar_recover (checking docstring)
"""

from typing import NamedTuple, List, Dict
from itertools import zip_longest


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
    
    valid_nucleotides = set('ACGTN')
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


def cigar_recover(sequence, cigar):
    """
    Recover the original sequence by interpreting the CIGAR string.
    
    CIGAR operations are handled as follows:
    - M: Keep the sequence as is.
    - D: Add '-' for deleted bases.
    - I: Convert bases to lowercase.
    - S: Keep the sequence as is (soft-clipped bases are retained).
    - H: Add '×' for hard-clipped bases.
    
    Args:
        sequence (str): The input sequence (e.g., 'CCTAGTCCAAACTGGATCTCTGCTGTCCCTG').
        cigar (str): The CIGAR string (e.g., '224H31M').
    
    Returns:
        str: The recovered sequence with modifications based on CIGAR operations.
    
    Examples:
        >>> sequence = 'CCTAGTCCAAACTGG'
        >>> cigar = '5M'  # Keep 5 bases as is
        >>> cigar_recover(sequence, cigar)
        'CCTAG'
        
        >>> cigar = '3H5M'  # Add 3 '×' then keep 5 bases
        >>> cigar_recover(sequence, cigar)
        '×××CCTAG'
        
        >>> cigar = '5M3D'  # Keep 5 bases, add 3 '-'
        >>> cigar_recover(sequence, cigar)
        'CCTAG---'
        
        >>> cigar = '5M3I'  # Keep 5 bases, lowercase 3 bases
        >>> cigar_recover(sequence, cigar)
        'CCTAGtcc'
        
        >>> cigar = '5M3S'  # Keep 5 bases, keep 3 bases (soft-clipped)
        >>> cigar_recover(sequence, cigar)
        'CCTAGTCC'
        
        >>> cigar = '2H3M2D2I3S'  # Combined: 2 '×', 3 bases, 2 '-', 2 lowercase, 3 bases
        >>> cigar_recover(sequence, cigar)
        '××CCT--aaACT'
        
        >>> sequence = 'CCTAGTCCAAACTGGATCTCTGCTGTCCCTG'
        >>> cigar = '224H31M'
        >>> cigar_recover(sequence, cigar)
        '××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××CCTAGTCCAAACTGGATCTCTGCTGTCCCTG'
    """
    import re
    
    # Parse CIGAR string into operations (e.g., [('224', 'H'), ('31', 'M')])
    cigar_ops = re.findall(r'(\d+)([HMSID])', cigar)
    
    result = []
    seq_pos = 0  # Position in the input sequence
    
    for count, op in cigar_ops:
        count = int(count)
        if op == 'H':
            # Add 'count' number of '×' for hard-clipped bases
            result.append('×' * count)
        elif op == 'D':
            # Add 'count' number of '-' for deleted bases
            result.append('-' * count)
        elif op == 'M' or op == 'S':
            # Keep 'count' bases from the sequence as is
            result.append(sequence[seq_pos:seq_pos + count])
            seq_pos += count
        elif op == 'I':
            # Convert 'count' bases to lowercase
            result.append(sequence[seq_pos:seq_pos + count].lower())
            seq_pos += count
    
    return ''.join(result)