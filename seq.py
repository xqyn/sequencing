'''
project: sequencing
april 22 2025 - XQ - Leiden UMC
modifying func for dna sequence
''' 

from typing import NamedTuple, List, Dict


# compute_base_percent --------------------------------------------------
def compute_base_count(seq_list: list,
                       bases: str = 'ATCGN-',
                       percentage: bool = False) -> dict:
    
    """
    calculate the count/percentage of each base at each position 
    in a   sequences list.
    
    Args:
        seq_list (list): list of DNA sequences (strings of equal length).
        bases (str): string of valid base characters (default: "ATCGN-").
        percentage (bool): whether to return as count or percentage (default: False) 
    
    Returns:
        dict: dict mapping each base to a list of percentages for each position.
    """
    if not seq_list:
        return {base: [] for base in bases}
    
    # assume length of sequence in seq_list are equal
    length = len(seq_list[0])
    
    # init dict for base count for each position 
    base_dict = {base: [0] * length for base in bases}
    
    # count occurences of each base from bases at each position
    for seq in seq_list:
            for pos, base in enumerate(seq):
                if base in bases:
                    base_dict[base][pos] += 1
    
    if percentage:
        return {base: [round(count / len(seq_list),5) * 100 for count in counts] for base, counts in base_dict.items()}
    return base_dict

