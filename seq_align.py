"""
project: sequencing
XQ - Leiden UMC
sequence modification: alignment
update: 2025
    - juni 17: pairwise alignment 
"""

import os
from Bio import Align


def align_pairwise(ref: str,
                   seq: str,
                   mode: str = 'local',
                   match_score: float = 1,
                   mismatch_score: float = 0,
                   open_gap_score: float = -0.9,
                   extend_gap_score: float = -2,
                   return_best: bool = True,     # Return aligned sequence + score or just score
                   return_alignment: bool = True):     # Return aligned sequence + score or just score
    
    """
    alignment method for pairwise sequence.

    Parameters
    ----------
    ref: 
        referemce sequence.
    seq: 
        sequence to align.
    mode: 
        alignment mode ('local' or 'global').
    match_score: 
        score for a match.
    mismatch_score: 
        score for a mismatch.
    open_gap_score: 
        gap open penalty.
    extend_gap_score: 
        gap extension penalty.
    return_alignment: 
        If True, returns (aligned sequence, score), else just score.
    print_alignment: 
        If True, prints the best alignment.
    
    Returns
    -------
        If return_alignment=True: (best_alignment_str, best_score)
        Else: best_score (float)
    
    Examples
    --------
    bc1 = 
    """
    # checking mode
    if mode not in ['local', 'global']:
        raise ValueError("Mode must be 'local' or 'global'")
    
    # create aligner model
    aligner = Align.PairwiseAligner()
    
    # setting aligne param    
    aligner.mode = mode
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score
    
    # performing alignement
    aligner_all = aligner.align(ref, seq)
    
    if return_best:
        best_align = None
        best_score = 0
        for align in aligner_all:
            if align.score > best_score:
                best_align = align
        return best_align
    else:
        return aligner_all

