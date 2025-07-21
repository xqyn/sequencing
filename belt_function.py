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
    - juli 21: add extract_read_sequence
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


# chr2_trace_edit.py--------------------------------------------------


# # --------------------------------------------------
# def trace_barcode(bc_sp_dict: dict[str, pd.DataFrame],
#                   bc_seq: str,
#                   condition: str = 'day') -> pd.DataFrame:
#     """
#     Trace editing changes for a specific barcode across experimental conditions.
    
#     Parameters
#     ----------
#     bc_sp_dict : dict of str -> pd.DataFrame
#         Dictionary mapping condition names (e.g., 'day1', 'day2') to DataFrames 
#         containing barcode-spacer mappings. Each DataFrame should have an index 
#         of barcode sequences and a 'spacer' column with corresponding spacer sequences.
#     bc_seq : str
#         The barcode sequence to trace across conditions.
#     condition : str, optional
#         The label name to use for the condition column in the result. Default is 'day'.
    
#     Returns
#     -------
#     pd.DataFrame
#         A DataFrame containing spacer sequences and their alignment scores and 
#         Hamming distances to the reference spacer for the specified barcode, 
#         across all conditions where the barcode is present.
    
#     Notes
#     -----
#     - The reference spacer is taken from the first condition where the barcode is found.
#     - Assumes `align_pairwise` returns an object with `.score`.
#     - Assumes `sequence.hamming(a, b)` returns the Hamming distance between two sequences.
    
#     Example
#     -------
#     >>> trace_barcode(bc_sp_dict, "ACGT123")
#     """
    
#     # collect all matching rows across conditions
#     df_trace = pd.DataFrame()
#     for cond in sorted(bc_sp_dict.keys()):
#         df_cond = bc_sp_dict[cond]
#         if bc_seq in df_cond.index:
#             df_cond = bc_sp_dict[cond][bc_sp_dict[cond].index == bc_seq]
#             df_cond.insert(2, condition, cond)
#             df_trace = pd.concat([df_trace, df_cond])
    
#     if df_trace.empty:
#         return df_trace  # Return empty DataFrame if barcode not found
    
#     # get original spacer and compute scores/distances
#     sp_ori = df_trace['spacer'].iloc[0]
#     df_trace['score'] = [seq_align.align_pairwise(sp_ori, spacer).score for spacer in df_trace['spacer']]
#     df_trace['hamm_dist'] = [sequence.hamming(sp_ori, spacer) for spacer in df_trace['spacer']]
    
#     return df_trace


# chr3_nam_check --------------------------------------------------
from Bio import SeqIO
import gzip
import sys


def extract_read_sequence(read1_file, read2_file, read_id, phred=False):
    """
    Extract sequences and optionally Phred scores for a given read ID from read1 and read2 FASTQ.gz files.
    
    Parameters:
    read1_file (str): Path to read1 FASTQ.gz file
    read2_file (str): Path to read2 FASTQ.gz file
    read_id (str): Read ID to search for (without '@')
    phred (bool): If True, return Phred quality scores; if False, return only sequences (default: False)
    
    Returns:
    tuple: (read1_sequence, read2_sequence) if phred=False, or
           (read1_sequence, read2_sequence, read1_phred, read2_phred) if phred=True.
           Returns (None, None) or (None, None, None, None) if not found or error occurs.
    """
    # Ensure read_id starts with '@' for FASTQ format
    read_id = read_id if read_id.startswith('@') else '@' + read_id
    
    read1_sequence = None
    read2_sequence = None
    read1_phred = None
    read2_phred = None
    
    try:
        # Parse read1 FASTQ.gz file
        with gzip.open(read1_file, 'rt') as handle1:  # 'rt' mode for text reading
            for record in SeqIO.parse(handle1, 'fastq'):
                if record.id == read_id.lstrip('@'):
                    read1_sequence = str(record.seq)
                    if phred:
                        read1_phred = record.letter_annotations["phred_quality"]
                    break
    except FileNotFoundError:
        print(f"Error: File {read1_file} not found.")
        return (None, None, None, None) if phred else (None, None)
    except Exception as e:
        print(f"Error reading {read1_file}: {str(e)}")
        return (None, None, None, None) if phred else (None, None)
    
    try:
        # Parse read2 FASTQ.gz file
        with gzip.open(read2_file, 'rt') as handle2:  # 'rt' mode for text reading
            for record in SeqIO.parse(handle2, 'fastq'):
                if record.id == read_id.lstrip('@'):
                    read2_sequence = str(record.seq)
                    if phred:
                        read2_phred = record.letter_annotations["phred_quality"]
                    break
    except FileNotFoundError:
        print(f"Error: File {read2_file} not found.")
        return (None, None, None, None) if phred else (None, None)
    except Exception as e:
        print(f"Error reading {read2_file}: {str(e)}")
        return (None, None, None, None) if phred else (None, None)
    
    return (read1_sequence, read2_sequence, read1_phred, read2_phred, record.id) if phred else (read1_sequence, read2_sequence)



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