'''
project: belt
XQ - Leiden UMC
Import global function
update:
    - aug 5 2024
    2025
    - april 
    - may 1
    - juni 17: add trace_barcode
    - juni 17: add hamming into collect_barcode
    - juli 21: add extract_seq_fastq
    - juli 23: add extract_pair_fastq, extract_bam, process_query/ries
    - juli 24-25: 
            add class fastq and ReadAlign
            copy `align_pairwise` from trace_barcode.py
            write out `process_pair_reads` for generate consensus read
    - juli 28: add local align for sp
    - juli 30: add barcode search
    - juli 31: add generate_consensus, decide on better base and score
            update adjust_sequence for position zero
note: combine trace_barcode and align_motif
''' 

# from Bio import Align
# from collections import Counter
# import numpy as np
# from itertools import zip_longest
# import sequence as sequence
# import seq_align as seq_align
# import re
# import belt_viz
# from Bio import SeqIO
# import gzip

import re
import gzip
import numpy as np
import pandas as pd
import pysam
from collections import Counter
from typing import List, Tuple, Optional, NamedTuple, Union
from itertools import zip_longest
from Bio import SeqIO
from Bio import Align
from Bio.Align import PairwiseAligner
import multiprocessing as mp
from functools import partial
import logging
import belt_viz

# Configure logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

# dna general --------------------------------------------------
def ascii_to_phred(quality_string, offset=33):
    """Convert an ASCII quality to Phred scores"""
    return [ord(char) - offset for char in quality_string]


def complement(seq: str, reverse: bool = False) -> str:
    """Return the complement (or reverse complement) of a DNA sequence."""
    mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
               'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}
    complemented = ''.join(mapping.get(base, base) for base in seq)
    return complemented[::-1] if reverse else complemented


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


# chr1_pre3_overview.py --------------------------------------------------
def count_uniq(df, column, return_list=False):
    """
    Returns either the number or list of unique values in `column` of `df` (those with count == 1).
    
    Args:
        df (pd.DataFrame): Input DataFrame
        column (str): Column name to analyze (index will be reset)
        return_list (bool): If True, return list of unique values instead of count
        
    Returns:
        int or list: Count or list of unique values
    """
    values = df.reset_index()[column]
    counts = Counter(values)
    
    uniques = [val for val, count in counts.items() if count == 1]
    
    return uniques if return_list else len(uniques)


# chr3_bam_check --------------------------------------------------
# ---BAM
def cigar_recover(sequence: str, cigar: str) -> str:
    """
    Recover original sequence by interpreting CIGAR string.
    
    CIGAR operations:
    - M: Keep sequence as is
    - D: Add '-' for deleted bases
    - I: Convert bases to lowercase
    - S: Keep sequence as is (soft-clipped)
    - H: Add '×' for hard-clipped bases
    
    Args:
        sequence: Input DNA sequence
        cigar: CIGAR string (e.g., '224H31M')
        
    Returns:
        Recovered sequence with CIGAR modifications applied
    Examples:
         >>> cigar = '2H3M2D2I3S'  # Combined: 2 '×', 3 bases, 2 '-', 2 lowercase, 3 bases
        >>> cigar_recover(sequence, cigar)
        '××CCT--aaACT'
    """
    if cigar is None:
        return sequence
    
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


def extract_bam(query_id: str, 
                bam_path: str) -> Tuple[List[str], List[str], List[str]]:
    """
    Extract sequences, CIGAR strings, and formatted names from a BAM file for a given query ID.
    
    Args:
        bam_path (str): Path to the BAM file.
        query_id (str): Query ID to filter reads.
    
    Returns:
        Tuple[List[str], List[str], List[str]]: Lists of sequences, CIGAR strings, and formatted names.
    """
    #logger.info(f"Processing BAM ID: {query_id}")
    seq_list, cigar_list, name_list = [], [], []
    
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for query in bam:
                if query.query_name == query_id:
                    seq_list.append(query.seq)
                    cigar_list.append(query.cigarstring)
                    name = 'read1_' if query.is_read1 else 'read2_'
                    reference_name = 'None' if query.reference_name is None else query.reference_name
                    cigarstring = 'None' if query.cigarstring is None else query.cigarstring
                    name_list.append(name + reference_name + '_' + cigarstring)
        
        # apply cigar_recover to each sequence-CIGAR pair
        seq_list = [cigar_recover(seq, cig)
                    for seq, cig in zip(seq_list, cigar_list)]
        return seq_list, cigar_list, name_list
    except Exception as e:
        logger.error(f"Error processing BAM: {query_id}: {str(e)}")
        return [], [], []


def reorder_reads(data: List[str], 
                  seq_list: List[str]) -> Tuple[List[str], List[str]]:
    """
    Reorder BAM read names and sequences based on predefined order.

    Args:
        data: List of read names.
        seq_list: List of corresponding sequences.

    Returns:
        Tuple of (reordered names, reordered sequences).
    """
    # Check if all elements in data start with 'ref', 'read1_', or 'read2_' and have valid suffixes
    valid_suffixes = ['start', 'inbetween', 'end', 'None']
    all_valid = all(
        x.startswith(('ref', 'read1_', 'read2_')) and 
        (x == 'ref' or x.split('_')[1] in valid_suffixes)
        for x in data
    )
    
    if not all_valid:
        return data, seq_list
    
    ref = [x for x in data if x == 'ref']
    read1 = sorted([x for x in data if x.startswith('read1_')], 
                  key=lambda x: valid_suffixes.index(x.split('_')[1]))
    read2 = sorted([x for x in data if x.startswith('read2_')], 
                  key=lambda x: valid_suffixes.index(x.split('_')[1]))
    reordered_data = ref + read1 + read2
    
    # Get indices of reordered data in original data
    indices = [data.index(x) for x in reordered_data]
    # Reorder seq_list using the same indices
    reordered_seq = [seq_list[i] for i in indices]
    
    return reordered_data, reordered_seq

# ---FastQC
class FastqRecord(NamedTuple):
    """Container for paired-end FASTQ record data."""
    read_id: str            # Read ID DNA 
    r1_seq: str             # Read 1 DNA sequence
    r1_phred: List[int]     # Read 1 Phred quality scores
    r2_seq: str             # Read 2 DNA sequence
    r2_phred: List[int]     # Read 2 Phred quality scores
    r2_seq_rev: str         # Reverse complement of Read 2 sequence
    r2_phred_rev: List[int] # Reversed Read 2 Phred scores
    

def extract_pair_fastq(read_id: str, 
                       read1_file: str, 
                       read2_file: str) -> Optional[FastqRecord]:
    """
    Extract paired-end FASTQ records for a given read ID.

    Args:
        read_id (str): Identifier of the sequence (with or without '@').
        read1_file (str): Path to the Read 1 FASTQ file (gzipped or plain).
        read2_file (str): Path to the Read 2 FASTQ file (gzipped or plain).

    Returns:
        Record or None if read not found or error occurs.
    """
    #logger.info(f"Processing FASTQ ID: {read_id}")
    # Ensure read_id starts with '@' for FASTQ format
    read_id = read_id if read_id.startswith('@') else '@' + read_id
    
    r1_id, r1_seq, r1_phred, r2_seq, r2_phred = None, None, None, None, None
    
    read_id = read_id if read_id.startswith('@') else f'@{read_id}'
    target_id = read_id.lstrip('@')
    
    def read_fastq_file(filepath: str) -> Tuple[Optional[str], Optional[List[int]]]:
        """Read a single FASTQ file and extract sequence and quality."""
        opener = gzip.open if filepath.endswith('.gz') else open
        
        try:
            with opener(filepath, 'rt') as file:
                for record in SeqIO.parse(file, 'fastq'):
                    if record.id == target_id:
                        return str(record.seq), record.letter_annotations["phred_quality"]
        except FileNotFoundError:
            logger.error(f"Error: File {filepath}: not found")
        except Exception as e:
            logger.error(f"Error reading {filepath}: {str(e)}")
            
        return None, None
    
    # Extract R1 and R2 data
    r1_seq, r1_phred = read_fastq_file(read1_file)
    r2_seq, r2_phred = read_fastq_file(read2_file)
    
    if not (r1_seq and r2_seq):
        print(f"Error: Read {read_id} not found in one or both files.")
        return None
    
    # Generate reverse complement and reversed quality scores
    r2_seq_rev = complement(r2_seq, reverse=True)
    r2_phred_rev = r2_phred[::-1]
    
    return FastqRecord(
        read_id=target_id,
        r1_seq=r1_seq,
        r1_phred=r1_phred,
        r2_seq=r2_seq,
        r2_phred=r2_phred,
        r2_seq_rev=r2_seq_rev,
        r2_phred_rev=r2_phred_rev
    )


# --- Alignment
class ReadAlign:
    def __init__(self, 
                 r1_seq_align: str, 
                 r1_phred_align: List[int], 
                 r2_seq_align: str,
                 r2_phred_align: List[int],
                 css_seq: str,
                 css_phred: List[int]):
        self.r1_seq_align = r1_seq_align
        self.r1_phred_align = r1_phred_align
        self.r2_seq_align = r2_seq_align
        self.r2_phred_align = r2_phred_align
        self.css_seq = css_seq
        self.css_phred = css_phred


def align_pairwise(ref: str,
                   seq: str,
                   mode: str = 'global',
                   match_score: float = 1,
                   mismatch_score: float = 0,
                   open_gap_score: float = -3,
                   extend_gap_score: float = -3,
                   return_best: bool = True):     # Return aligned sequence + score or just score
    
    """
    alignment method for pairwise sequence.

    args:
        ref: referemce sequence.
        seq: sequence to align.
        mode: alignment mode ('local' or 'global').
        match_score: score for a match.
        mismatch_score: score for a mismatch.
        open_gap_score: gap open penalty.
        extend_gap_score: gap extension penalty.
        return_alignment: If True, returns (aligned sequence, score), else just score.
        print_alignment: If True, prints the best alignment.
    
    returns
        if return_alignment=True: (best_alignment_str, best_score)
        else: best_score (float)

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
    #aligner.target_left_gap_score = 0  
    aligner.target_left_gap_score = -5  # High penalty to prevent gaps at target start
    aligner.target_right_gap_score = -30
    aligner.query_left_open_gap_score = 0  # Allow gap at query start
    
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


# def adjust_sequence(ref: str, 
#                     seq: str, 
#                     seq_phred: List[int], 
#                     ref_start: int = 0, 
#                     seq_start: int = 0, 
#                     shift: Optional[int] = None, 
#                     compensate: bool = False) -> Tuple[str, List[int]]:
#     """
#     Adjusts a sequence and its Phred scores by shifting and/or compensating length.
    
#     Args:
#         ref (str): Reference sequence.
#         seq (str): Input sequence to adjust.
#         seq_phred (list): Phred scores corresponding to seq.
#         ref_start (int, optional): Reference sequence start index. Defaults to 0.
#         seq_start (int, optional): Input sequence start index. Defaults to 0.
#         compensate (bool, optional): If True, adjusts seq length to match ref. Defaults to False.
    
#     Returns:
#         tuple: Adjusted (sequence, Phred scores).
    
#     Raises:
#         ValueError: If seq and seq_phred lengths do not match.
#     """
#     if len(seq) != len(seq_phred):
#         raise ValueError("Sequence and Phred scores must have equal length")
    
#     # Use provided shift or calculate from seq_start - ref_start
#     shift = shift if shift is not None else seq_start - ref_start
    
#     # apply shift
#     if shift != 0:
#         if shift > 0:
#             seq = seq[shift:]
#             seq_phred = seq_phred[shift:]
#         else:
#             seq = ref[:abs(shift)].lower() + seq
#             seq_phred = [0] * abs(shift) + seq_phred
    
#     # Apply compensation if requested
#     if compensate:
#         len_diff = len(seq) - len(ref)
#         if len_diff > 0:
#             seq = seq[:len(ref)]
#             seq_phred = seq_phred[:len(ref)]
#         elif len_diff < 0:
#             seq += ref[-abs(len_diff):].lower()
#             seq_phred += [0] * abs(len_diff)
    
#     return seq, seq_phred


def adjust_sequence(ref: str, 
                    seq: str, 
                    seq_phred: List[int], 
                    pos_motif_ref: int = 0,
                    pos_align_seq: int = 0, 
                    shift: Optional[int] = None, 
                    compensate: bool = False) -> Tuple[str, List[int]]:
    """
    Adjusts a sequence and its Phred scores by shifting and/or compensating length.
    
    Args:
        ref (str): Reference sequence.
        seq (str): Input sequence to adjust.
        seq_phred (list): Phred scores corresponding to seq.
        ref_start (int, optional): Reference sequence start index. Defaults to 0.
        seq_start (int, optional): Input sequence start index. Defaults to 0.
        compensate (bool, optional): If True, adjusts seq length to match ref. Defaults to False.
    
    Returns:
        tuple: Adjusted (sequence, Phred scores).
    
    Raises:
        ValueError: If seq and seq_phred lengths do not match.
    Describe:
        [pos]: the position of start
        shift: the different in two seq
        ×××××: motif
        [=====]: barcode/spacer of target
        
        pos_motif_ref = 0
        shift > 0
            ref: |------------------------------------------|
            seq: |shift-------------------------------------|
        shift < 0
            ref: |     -------------------------------------|
            seq: |shift-------------------------------------|
        pos_motif_ref > 0 (abbrx. pos)
        shift > 0
            ref: |----------[pos]×××××[=====]---------------|
            seq: |----------shift[pos]×××××[=====]----------|
        shift < 0
            ref: |----------[pos]×××××[=====]---------------|
            seq: |-----shift[pos]×××××[=====]---------------|
    """
    if len(seq) != len(seq_phred):
        raise ValueError("Sequence and Phred scores must have equal length")
    
    # Use provided shift or calculate from seq_start - ref_start
    shift = shift if shift is not None else pos_align_seq - pos_motif_ref
    pos_align_seq = pos_motif_ref+shift
    
    # apply shift
    if shift != 0:
        if shift > 0:
            seq = seq[:pos_motif_ref] + seq[pos_align_seq:]
            seq_phred = seq_phred[:pos_motif_ref] + seq_phred[pos_align_seq:]
        else:
            if pos_motif_ref != 0:
                seq = seq[:pos_align_seq] + ref[pos_align_seq:pos_motif_ref].lower() + seq[pos_align_seq:]
                seq_phred = seq_phred[:pos_align_seq] + [0]*abs(shift) + seq_phred[pos_align_seq:]
            else:
                seq = ref[:abs(shift)].lower() + seq
                seq_phred = [0]*abs(shift) + seq_phred
                          
    # Apply compensation if requested
    if compensate:
        len_diff = len(seq) - len(ref)
        if len_diff > 0:
            seq = seq[:len(ref)]
            seq_phred = seq_phred[:len(ref)]
        elif len_diff < 0:
            seq += ref[-abs(len_diff):].lower()
            seq_phred += [0] * abs(len_diff)
    
    return seq, seq_phred


def generate_consensus(r1_seq_align: str,
                       r2_seq_align: str, 
                       r1_phred_align: str, 
                       r2_phred_align: str, 
                       ref: str) -> tuple[str, str]:
     
    middle = round(len(r1_seq_align)/2) - 1
    phred_r1_1st = sum(r1_phred_align[:middle]) / len(r1_phred_align[:middle])
    phred_r2_1st = sum(r2_phred_align[:middle]) / len(r2_phred_align[:middle])
    phred_r1_2nd = sum(r1_phred_align[middle:]) / len(r1_phred_align[middle:])
    phred_r2_2nd = sum(r2_phred_align[middle:]) / len(r2_phred_align[middle:])
    
    css_seq, css_phred = [], []
    for ind, (r1_score, r2_score) in enumerate(zip(r1_phred_align, r2_phred_align)):
        r1_base, r2_base = r1_seq_align[ind], r2_seq_align[ind]            
        if r1_score != r2_score:
            base, score = (r1_base, r1_score) if r1_score > r2_score else (r2_base, r2_score)
        else:
            # Equal scores — resolve by base match or context
            base = (
                r1_base if r1_base == r2_base or r1_base == ref[ind] else
                r2_base if r2_base == ref[ind] else
                (
                    r1_base if (ind < middle and phred_r1_1st > phred_r2_1st) or 
                    (ind > middle and phred_r1_2nd > phred_r2_2nd) else
                    r2_base if (ind < middle and phred_r1_1st < phred_r2_1st) or 
                    (ind > middle and phred_r1_2nd < phred_r2_2nd) else
                    'Q')
            )
            score = r1_score if base == r1_base else r2_score if base == r2_base else 0
        
        css_seq.append(base)
        css_phred.append(score)
    
    return ''.join(css_seq), css_phred

# css_seq, css_phred = generate_consensus(r1_seq_align,
#                                                 r2_seq_align,
#                                                 r1_phred_align,
#                                                 r2_phred_align,
#                                                 ref)

def process_pair_reads(ref: str, 
                       read_pairs: List[Tuple[str, List[int]]], 
                       query_id: str = 'test') -> Optional[ReadAlign]:
    try:
        r1_seq_align, r1_phred_align, r2_seq_align, r2_phred_align = None, None, None, None
        # Process each read pair
        for pair_idx, (seq, seq_phred) in enumerate(read_pairs, 1):
            #print(f"Processing read pair {pair_idx}...")
            
            # step 1: initial alignment and check for sequencing shift
            aligner = align_pairwise(ref, seq, mode='global')
            start_ref = aligner.indices[0][0]
            start_seq = aligner.indices[1][0]
            shift = start_seq - start_ref
            
            if shift != 0:
                seq, seq_phred = adjust_sequence(ref, 
                                                seq, 
                                                seq_phred, 
                                                shift=shift, 
                                                compensate=False)
                aligner = align_pairwise(ref, seq, mode='global')
            else:
                # no shift detected, proceeding with original sequence
                pass
            
            # step 2: check for insertions and remove them
            ref_align = aligner.indices[0]
            insert = np.where(ref_align < 0)[0].tolist()
            if insert:
                seq = ''.join([c for i, c in enumerate(seq) if i not in insert])
                seq_phred = [v for i, v in enumerate(seq_phred) if i not in insert]
            else:
                pass
            
            # step 3: final alignment and compute aligned Phred scores
            aligner = align_pairwise(ref, seq, mode='global')
            seq = aligner[1]
            seq_phred = [0 if pos < 0 else seq_phred[pos] for pos in aligner.indices[1]]
            
            # step 4: add or trim artificial bases if len difference
            seq, seq_phred = adjust_sequence(ref,
                                            seq, 
                                            seq_phred, 
                                            shift=0, 
                                            compensate=True)
                
            # store results for this read pair
            if pair_idx == 1:
                r1_seq_align, r1_phred_align = seq, seq_phred
                if len(r1_seq_align) != len(r1_phred_align):
                    raise ValueError(f"Length mismatch for read 1: sequence length {len(r1_seq_align)} != PHRED length {len(r1_phred_align)}")
            else:
                r2_seq_align, r2_phred_align = seq, seq_phred
                if len(r2_seq_align) != len(r2_phred_align):
                    raise ValueError(f"Length mismatch for read 2: sequence length {len(r2_seq_align)} != PHRED length {len(r2_phred_align)}")
        
        
        # step 5: generate consensus sequence and Phred scores
        logger.info("Generating consensus sequence...")
        css_seq, css_phred = generate_consensus(r1_seq_align,
                                                r2_seq_align,
                                                r1_phred_align,
                                                r2_phred_align,
                                                ref)
        
        return ReadAlign(r1_seq_align, 
                        r1_phred_align, 
                        r2_seq_align, 
                        r2_phred_align, 
                        css_seq, 
                        css_phred)
    except Exception as e:
            logger.error(f"Error processing query {query_id}: {str(e)}")
            return  # Continue to next query



# --- process query
def process_query(index: int,
                  query_id: str,
                  read1_path: str,
                  read2_path: str,
                  read1_full_path: str,
                  read2_full_path: str,
                  bam_path: str,
                  ref: str,
                  motif_sp_down: str = 'CGGGGTTAGA',
                  motif_sp_down_pos: int = 100,
                  motif_bp_up: str  = 'GCTTTTTTTT',
                  motif_bp_up_pos: int = 176,
                  shift_base: int = 3,
                  sample: str = 'test',
                  fig_dir: str = './figure/test/'
                  ) -> None:
    """
    Analyzes FASTQ and BAM sequences, generating visualization grids with optional multiprocessing.
    
    Args:
        index: Query index for logging.
        query_id: Query identifier.
        read1_path: Path to read1 FASTQ file.
        read2_path: Path to read2 FASTQ file.
        read1_full_path: Path to full read1 FASTQ file.
        read2_full_path: Path to full read2 FASTQ file.
        bam_path: Path to BAM file.
        ref: Reference sequence.
        motif: Motif sequence for alignment check.
        motif_sp_down_pos: Expected position of motif.
        shift_base: Maximum allowed shift before correction.
        sample: Sample name for output files.
        fig_dir: Directory for output figures.
    
    Steps:
        1. Iterate through indices (in parallel if use_multiprocessing=True).
        2. Clear global variables (read1, r2_seq, r1_phred, r2_phred, id).
        3. Extract FASTQ sequences and PHRED scores for given query ID.
        4. Compute reverse complement of r2_seq and reverse its PHRED scores.
        5. Plot FASTQ grid with reference, reads, and PHRED scores.
        6. Process BAM data: collect reads, sequences, CIGAR strings, and names.
        7. Recover sequences using CIGAR strings and insert reference.
        8. Extract numbers from CIGAR strings and combine with default x-labels.
        9. Reorder reads and plot BAM mapping grid.
    """
    # remove variables
    for var in ['read1', 'r2_seq', 'r1_phred', 'r2_phred', 'read_id']:
        try:
            del globals()[var]
        except KeyError:
            pass
    
    logger.info(f"Processing query {index}: {query_id}")
    
    # --- BAM: check mapping
    x_labels = [0, 1, 20, 73, 84, 99, 100, 120, 130, 140, 150, 160, 170, 177, 186, 223, 224, 250, 255]
    bam_seq_list, bam_name_list, numbers = [], [], []
    logger.info(f"Processing BAM file...")
    bam_seq_list, bam_cigar_list, bam_name_list = extract_bam(query_id, bam_path)
    # re-order start-inbetween-end
    bam_name_list, bam_seq_list = reorder_reads(bam_name_list, bam_seq_list)
    # find the number in bam for easier to follow and add xlabels
    numbers = [int(num) for item in bam_cigar_list if item for num in re.findall(r'\d+', item)]
    
    # for item in bam_cigar_list:
    #     if item is not None:  # Skip None values
    #         found_numbers = re.findall(r'\d+', item)
    #         numbers.extend(int(num) for num in found_numbers)
    x_labels = sorted(list(set(numbers + x_labels)))
    
    # --- FASTQC: check raw read
    logger.info(f"Processing FASTQC file...")
    fq = extract_pair_fastq(query_id, read1_path, read2_path)
    if not fq:
        logger.error(f"Failed to extract FASTQ for query {query_id}")
        return
    
    r1_phred, r1_seq = fq.r1_phred, fq.r1_seq
    r2_seq_rev, r2_phred_rev = fq.r2_seq_rev, fq.r2_phred_rev
    
    # --- ALIGNMENTS:
    # -- spacer
    # check if the read have hugh shift
    # read 1 
    r1_seq_sp, r1_phred_sp = r1_seq, r1_phred
    motif_align = align_pairwise(r1_seq_sp, motif_sp_down, mode='local')
    shift_r1 = motif_align.indices[0][0] - motif_sp_down_pos
    
    if abs(shift_r1) > 80:
        logger.warning(f"Huge shift {shift_r1} in Read 1, skipping index {index}")
        return
    if abs(shift_r1) > shift_base:
        logger.warning(f"Correcting shift {shift_r1} in Read 1")
        fq_full = extract_pair_fastq(query_id, read1_full_path, read2_full_path)
        logger.info(f"Processing FASTQC finalizing...")
        if fq_full:
            logger.info(f"Adjusting new sequence...")
            r1_seq_sp, r1_phred_sp = adjust_sequence(ref, 
                                                     fq_full.r1_seq, 
                                                     fq_full.r1_phred, 
                                                     shift=shift_r1, 
                                                     compensate=True)
    
    # read 2 
    r2_seq_rev_sp, r2_phred_rev_sp = r2_seq_rev, r2_phred_rev
    motif_align = align_pairwise(r2_seq_rev_sp, motif_sp_down, mode='local')
    shift_r2 = motif_align.indices[0][0] - motif_sp_down_pos
    
    if abs(shift_r2) > 80:
        logger.warning(f"Huge shift {shift_r2} in Read 2, skipping index {index}")
        return
    if abs(shift_r2) > shift_base:
        logger.warning(f"Correcting shift {shift_r2} in Read 2")
        fq_full = extract_pair_fastq(query_id, read1_full_path, read2_full_path)
        logger.info(f"Processing FASTQC finalizing...")
        if fq_full:
            logger.info(f"Adjusting new sequence...")
            r2_seq_rev_sp, r2_phred_rev_sp = adjust_sequence(ref, 
                                                             fq_full.r2_seq_rev, 
                                                             fq_full.r2_phred_rev, 
                                                             shift=shift_r2 + 38, 
                                                             compensate=True)
    
    # --- barcode 
    # read 1
    motif_align = align_pairwise(r1_seq_sp, motif_bp_up, mode='local')
    r1_seq_bc, r1_phred_bc = adjust_sequence(ref=ref,
                                             seq=r1_seq_sp, 
                                             seq_phred=r1_phred_sp,
                                             pos_motif_ref=motif_bp_up_pos,
                                             pos_align_seq=motif_align.indices[0][0],
                                             compensate=True)
    
    # read 2
    motif_align = align_pairwise(r2_seq_rev_sp, motif_bp_up, mode='local')
    r2_seq_rev_bc, r2_phred_rev_bc = adjust_sequence(ref=ref,
                                             seq=r2_seq_rev_sp, 
                                             seq_phred=r2_phred_rev_sp,
                                             pos_motif_ref=motif_bp_up_pos,
                                             pos_align_seq=motif_align.indices[0][0],
                                             compensate=True)

    # ---final
    logger.info(f"Processing pair reads...")
    pair = process_pair_reads(ref, 
                              [(r1_seq_bc, r1_phred_bc), (r2_seq_rev_bc,r2_phred_rev_bc)], 
                              query_id)
    
    
    # --- plots
    logger.info(f"Plotting {index}...")
    section = '-'*len(ref)
    seq_list = bam_seq_list + [section, ref,
                               fq.r1_phred, fq.r1_seq, fq.r2_seq_rev, fq.r2_phred_rev,
                               section, ref,
                               r1_phred_sp, r1_seq_sp, r2_seq_rev_sp, r2_phred_rev_sp, 
                               section, ref, 
                               r1_phred_bc, r1_seq_bc, r2_seq_rev_bc, r2_phred_rev_bc,
                               section, ref,
                               pair.r1_phred_align, pair.r1_seq_align, pair.r2_seq_align, pair.r2_phred_align,
                               section, ref,
                               pair.css_seq, pair.css_phred]
    
    y_labels= bam_name_list + ['', 'ref', 
                              '255_phred_r1', '255_seq_r1', '255_rev_seq_r2', '255_phred_rev_r2',
                              '','ref', 
                              'r1_phred_sp', 'r1_seq_sp', 'r2_seq_rev_sp',  'r2_phred_rev_sp', 
                              '','ref', 
                              'r1_phred_bc', 'r1_seq_bc', 'r2_seq_rev_bc', 'r2_phred_rev_bc',
                              '', 'ref',
                              'align_phred_r1', 'align_r1', 'align_r2', 'aling_phred_r2',
                              '', 'ref',
                              'seq_consensus','phred_consensus']
    # --- plot
    try:
        belt_viz.plot_base_grid(
            seq_list,
            x_labels=x_labels,
            y_labels=y_labels,
            height=6,
            width=30,
            fig_dir=fig_dir,
            fig_name=f"{sample}_{index}_align_combine",
            plot_title=f"{sample}_{index}_align_combine",
            fsize=6,
            artificial=True
        )
    except Exception as e:
        logger.error(f"Error generating visualization for {query_id}: {str(e)}")


# ---multiprocessing
def process_queries(index_query_pairs,
                    read1_path,
                    read2_path, 
                    read1_full_path,
                    read2_full_path,
                    bam_path,
                    ref,
                    motif_sp_down = 'CGGGGTTAGA',
                    motif_sp_down_pos = 100,
                    motif_bp_up = 'GCTTTTTTTT',
                    motif_bp_up_pos = 176,
                    shift_base=10,
                    sample='test', 
                    fig_dir='./figure/test/',
                    num_processes=None):
    """
    Process multiple (index, query_id) pairs in parallel using multiprocessing.
    
    Args:
        index_query_pairs (list): List of tuples [(index1, query_id1), (index2, query_id2), ...]
        read1_path (str): Path to read1 FASTQ file
        read2_path (str): Path to read2 FASTQ file
        bam_path (str): Path to output BAM file
        ref (str): Path to reference genome
        sample (str): Sample name (default: 'test')
        fig_dir (str): Directory for output figures (default: './figure/test/')
        num_processes (int): Number of processes to use (default: None, uses CPU count)
    
    Returns:
        None
    """
    # Determine number of processes
    if num_processes is None:
        num_processes = mp.cpu_count()
    
    # Create a partial function with fixed arguments
    process_func = partial(
        process_query,
        read1_path=read1_path,
        read2_path=read2_path,
        read1_full_path=read1_full_path,
        read2_full_path=read2_full_path,
        bam_path=bam_path,
        ref=ref,
        motif_sp_down = motif_sp_down,
        motif_sp_down_pos = motif_sp_down_pos,
        motif_bp_up = motif_bp_up,
        motif_bp_up_pos = motif_bp_up_pos,
        shift_base=shift_base,
        sample=sample,
        fig_dir=fig_dir
    )
    
    # Create a process pool
    with mp.Pool(processes=num_processes) as pool:
        results = [pool.apply_async(process_func, args=(pair[0], pair[1])) for pair in index_query_pairs]
        for r in results:
            try:
                r.wait()  # Wait for each process to complete
            except Exception as e:
                print(f"Process failed: {e}")
                continue  # Continue with next process
        # # Map the process function to index_query_pairs
        # pool.starmap(process_func, [(pair[0], pair[1]) for pair in index_query_pairs])
    
    print("Processes completed.")
