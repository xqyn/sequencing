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
note: combine trace_barcode and align_motif
''' 

# import sys
# import yaml
# from box import Box
# # # import variables:
# cfg = Box(yaml.safe_load(open('config.yaml')), default_box=True)
# # import utils from path:
# if cfg.dir.utils not in sys.path:
#     sys.path.append(cfg.dir.utils)

from Bio import Align
from collections import Counter
import numpy as np
import pandas as pd

import sequence as sequence
import seq_align as seq_align
import re
import belt_viz
from Bio import SeqIO
import gzip
import sys


# dna general --------------------------------------------------

def ascii_to_phred(quality_string, offset=33):
    """Convert an ASCII quality to Phred scores"""
    return [ord(char) - offset for char in quality_string]


def rever_complement(seq):
    """Return the reverse complement of a sequence"""
    return complement(seq, reverse=True)


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

# ---FastQC
from typing import NamedTuple
from Bio import SeqIO
import gzip

class fastq(NamedTuple):
    """
    FASTQ record storage.

    Attributes:
        read_id (str): Sequence identifier.
        r1_seq (str): Read 1 DNA sequence.
        r1_phred (list): Read 1 Phred quality scores.
        r2_seq (str): Read 2 DNA sequence.
        r2_phred (list): Read 2 Phred quality scores.
        r2_seq_rev (str): Reverse complement of Read 2 sequence.
        r2_phred_rev (list): Reversed Read 2 Phred scores.
    """
    read_id: str
    r1_seq: str
    r1_phred: list
    r2_seq: str
    r2_phred: list
    r2_seq_rev: str
    r2_phred_rev: list
    

def extract_pair_fastq(read_id,
                      read1_file,
                      read2_file):
    """
    Extract paired-end FASTQ records for a given read ID.

    Args:
        read_id (str): Identifier of the sequence (with or without '@').
        read1_file (str): Path to the Read 1 FASTQ file (gzipped or plain).
        read2_file (str): Path to the Read 2 FASTQ file (gzipped or plain).

    Returns:
        fastq: NamedTuple containing read ID, R1 sequence, R1 Phred scores,
               R2 sequence, R2 Phred scores, R2 reverse complement sequence,
               and R2 reverse Phred scores. Returns None if read not found or error occurs.
    """
    # Ensure read_id starts with '@' for FASTQ format
    print(f"Processing fastq ID: {read_id}")
    read_id = read_id if read_id.startswith('@') else '@' + read_id
    
    r1_seq = r1_phred = r2_seq = r2_phred = None
    r1_id = None
    
    # Parse R1 FASTQ
    #print(f"Opening R1 file...")
    opener = gzip.open if read1_file.endswith('.gz') else open
    try:
        with opener(read1_file, 'rt') as queryR1:
            for r1_query in SeqIO.parse(queryR1, 'fastq'):
                if r1_query.id == read_id.lstrip('@'):
                    r1_id = r1_query.id
                    r1_seq = str(r1_query.seq)
                    r1_phred = r1_query.letter_annotations["phred_quality"]
                    break
    except FileNotFoundError:
        print(f"Error: File {read1_file} not found.")
        return None
    except Exception as e:
        print(f"Error reading {read1_file}: {str(e)}")
        return None
    
    #print(f"Opening R2 file...")
    # Parse R2 FASTQ
    opener = gzip.open if read2_file.endswith('.gz') else open
    try:
        with opener(read2_file, 'rt') as queryR2:
            for r2_query in SeqIO.parse(queryR2, 'fastq'):
                if r2_query.id == read_id.lstrip('@'):
                    r2_seq = str(r2_query.seq)
                    r2_phred = r2_query.letter_annotations["phred_quality"]
                    break
    except FileNotFoundError:
        print(f"Error: File {read2_file} not found.")
        return None
    except Exception as e:
        print(f"Error reading {read2_file}: {str(e)}")
        return None
    
    if not (r1_seq and r2_seq):
        print(f"Error: Read {read_id} not found in one or both files.")
        return None
    
    # reverse complement of r2_seq and its phredscore
    def complement(seq, reverse=False):
        """Return the complement of a sequence"""
        mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
                   'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}
        complemented = ''.join([mapping.get(base, base) for base in seq])
        return complemented[::-1] if reverse else complemented
    
    r2_seq_rev = complement(r2_seq, reverse=True)
    r2_phred_rev = r2_phred[::-1]
    
    return fastq(
        read_id=r1_id,
        r1_seq=r1_seq,
        r1_phred=r1_phred,
        r2_seq=r2_seq,
        r2_phred=r2_phred,
        r2_seq_rev=r2_seq_rev,
        r2_phred_rev=r2_phred_rev
    )

def align_pairwise(ref: str,
                   seq: str,
                   mode: str = 'local',
                   match_score: float = 1,
                   mismatch_score: float = 0,
                   open_gap_score: float = -0.9,
                   extend_gap_score: float = -2,
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
    aligner.target_left_gap_score = 0  
    aligner.target_left_gap_score = -10  # High penalty to prevent gaps at target start
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

# ---Making consensus
class ReadAlign:
    def __init__(self, 
                 r1_seq_align, 
                 r1_phred_align, 
                 r2_seq_align, 
                 r2_phred_align, 
                 css_seq, 
                 css_phred):
        self.r1_seq_align = r1_seq_align
        self.r1_phred_align = r1_phred_align
        self.r2_seq_align = r2_seq_align
        self.r2_phred_align = r2_phred_align
        self.css_seq = css_seq
        self.css_phred = css_phred

def process_pair_reads(ref, read_pairs):
    # Process each read pair
    for pair_idx, (seq, seq_phred) in enumerate(read_pairs, 1):
        #print(f"\nProcessing read pair {pair_idx}...")
        
        # step 1: initial alignment and check for sequencing shift
        aligner = align_pairwise(ref, seq, mode='local')
        start_ref = aligner.indices[0][0]
        start_seq = aligner.indices[1][0]
        shift = start_seq - start_ref
        # adjust sequence and Phred scores based on shift
        if shift > 0:
            # shift > 0: trimming first {shift} bases from sequence and Phred scores
            seq = seq[shift:]
            seq_phred = seq_phred[shift:]
            aligner = align_pairwise(ref, seq, mode='local')
        elif shift < 0:
            seq = ref[:abs(shift)] + seq
            seq_phred = [0] * abs(shift) + seq_phred
            aligner = align_pairwise(ref, seq, mode='local')
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
        aligner = align_pairwise(ref, seq, mode='local')
        seq = aligner[1]
        phred_align = [0 if pos < 0 else seq_phred[pos] for pos in aligner.indices[1]]
        
        # step 4: add artificial bases if short len
        if len(seq) < len(ref):
            add = len(ref) - len(seq)
            seq += ref[-add:]
            phred_align += [0] * add
            
            
        # store results for this read pair
        if pair_idx == 1:
            r1_seq_align, r1_phred_align = seq, phred_align
        else:
            r2_seq_align, r2_phred_align = seq, phred_align
    
    # step 5: generate consensus sequence and Phred scores
    print("Generating consensus sequence...")
    css_seq = []
    css_phred = []
    for ind, (r1_base_score, r2_base_score) in enumerate(zip(r1_phred_align, r2_phred_align)):
        mean = int((r1_base_score + r2_base_score) / 2)
        r1_base, r2_base = r1_seq_align[ind], r2_seq_align[ind]
        
        if r1_base_score > r2_base_score:
            base = r1_base
        elif r2_base_score > r1_base_score:
            base = r2_base
        else:
            base = (r1_base if r1_base == r2_base or r1_base == ref[ind] else
                    r2_base if r2_base == ref[ind] else 'Q')
        
        css_seq.append(base)
        css_phred.append(mean if base != 'Q' else 0)
    
    # Return results as a single object
    
    return ReadAlign(r1_seq_align, 
                     r1_phred_align, 
                     r2_seq_align, 
                     r2_phred_align, 
                     ''.join(css_seq), 
                     css_phred)

# ---BAM
import pysam
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


def extract_bam(
        query_id: str,
        bam_path: str):
    """
    Extract sequences, CIGAR strings, and formatted names from a BAM file for a given query ID.
    
    Args:
        bam_path (str): Path to the BAM file.
        query_id (str): Query ID to filter reads.
    
    Returns:
        Tuple[List[str], List[str], List[str]]: Lists of sequences, CIGAR strings, and formatted names.
    """
    print(f"Processing bam ID: {query_id}")
    seq_list = []
    cigar_list = []
    name_list = []
    
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for query in bam:
            if query.query_name == query_id:
                seq_list.append(query.seq)
                cigar_list.append(query.cigarstring)
                name = 'read1_' if query.is_read1 else 'read2_'
                reference_name = 'None' if query.reference_name is None else query.reference_name
                cigarstring = 'None' if query.cigarstring is None else query.cigarstring
                name_list.append(name + reference_name + '_' + cigarstring)
    
    # Apply cigar_recover to each sequence-CIGAR pair
    seq_list = [cigar_recover(seq, cig) for seq, cig in zip(seq_list, cigar_list)]
    
    return seq_list, cigar_list, name_list


# def reorder_reads(data, seq_list):
#     ref = [x for x in data if x == 'ref']
#     read1 = sorted([x for x in data if x.startswith('read1_')], key=lambda x: ['start', 'inbetween', 'end', 'None'].index(x.split('_')[1]))
#     r2_seq = sorted([x for x in data if x.startswith('read2_')], key=lambda x: ['start', 'inbetween', 'end', 'None'].index(x.split('_')[1]))
#     reordered_data = ref + read1 + r2_seq
    
#     # Get indices of reordered data in original data
#     indices = [data.index(x) for x in reordered_data]
#     # Reorder seq_list using the same indices
#     reordered_seq = [seq_list[i] for i in indices]
    
#     return reordered_data, reordered_seq

def reorder_reads(data, seq_list):
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

# --- process query
def process_query(index, 
                  query_id,
                  read1_path,
                  read2_path, 
                  bam_path,
                  ref, 
                  sample='test', 
                  fig_dir='./figure/test/'):
    """
    Analyzes FASTQ and BAM sequences, generating visualization grids with optional multiprocessing.
    
    Args:
        indices (list): List of indices to process.
        query_id_list (list): List of query IDs.
        read1_path (str): Path to read1 FASTQ file.
        r2_seq_path (str): Path to r2_seq FASTQ file.
        query_list (list): List of query objects for BAM analysis.
        ref (str): Reference sequence.
        sample (str): Sample name for file naming.
        fig_dir (str): Directory to save output figures.
    
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
    
    # --- BAM: check mapping
    print(f'Run BAM: {index}')
    bam_seq_list = []
    bam_name_list = []
    x_labels = [0, 1, 20, 73, 84, 99, 186, 223, 250, 255]
    bam_seq_list, bam_cigar_list, bam_name_list = extract_bam(query_id,
                                                  bam_path)
    
    # re-order start-inbetween-end
    bam_name_list, bam_seq_list = reorder_reads(bam_name_list, bam_seq_list)
    numbers = []
    for item in bam_cigar_list:
        if item is not None:  # Skip None values
            found_numbers = re.findall(r'\d+', item)
            numbers.extend(int(num) for num in found_numbers)
    
    x_labels = sorted(list(set(numbers + x_labels)))
    
    # --- FASTQC: check raw read
    print(f"Run FASTQ: {index}")
    fq = extract_pair_fastq(query_id, read1_path, read2_path)

    r1_phred = fq.r1_phred
    r1_seq = fq.r1_seq
    r2_seq_rev = fq.r2_seq_rev
    r2_phred_rev = fq.r2_phred_rev
    
    # --- ALIGNMENTS:
    pair = process_pair_reads(ref, [(r1_seq, r1_phred),(r2_seq_rev,r2_phred_rev)])


    # --- plots
    section = '-'*len(ref)
    seq_list = bam_seq_list + [section, ref,
                               r1_phred, r1_seq, r2_seq_rev, r2_phred_rev,
                               section,
                               pair.r1_phred_align, pair.r1_seq_align, pair.r2_seq_align, pair.r2_phred_align,
                               section,
                               ref,
                               pair.css_seq,
                               pair.css_phred]

    ylabels= bam_name_list + ['',
             'ref', 'phred_r1', 'seq_r1', 'rev_seq_r2', 'phred_rev_r2',
             '',
             'align_phred_r1', 'align_r1', 'align_r2', 'aling_phred_r2',
             '',
             'ref',
             'seq_consensus','phred_consensus']
    
    # --- plot
    belt_viz.plot_base_grid(
        seq_list,
        x_labels=x_labels,
        y_labels=ylabels,
        height=3,
        width=30,
        fig_dir=fig_dir,
        fig_name=f"{sample}_{index}_align_combine",
        plot_title=f"{sample}_{index}_align_combine",
        fsize=6,
        artificial=True
    )


# ---multiprocessing
import multiprocessing as mp
from functools import partial


def process_queries(index_query_pairs,
                    read1_path,
                    read2_path, 
                    bam_path,
                    ref, 
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
        bam_path=bam_path,
        ref=ref,
        sample=sample,
        fig_dir=fig_dir
    )
    
    # Create a process pool
    with mp.Pool(processes=num_processes) as pool:
        # Map the process function to index_query_pairs
        pool.starmap(process_func, [(pair[0], pair[1]) for pair in index_query_pairs])
    
    print("Processes completed.")
