'''
project: belt
XQ - Leiden UMC
Import global function
update:
    - aug 5 2024
    - april 2025
    - may 1 2025
    - juni 17 25: add trace_barcode
    - juni 17 25: add hamming into collect_barcode
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