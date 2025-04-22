#!/bin/bash python

'''
project: sequencing
march 14 2025 - XQ - Leiden UMC
generating barcode variants
''' 


print('Importing barcode_utils.py')

import pandas as pd
from collections import Counter
from itertools import combinations, product


def generate_barcode_variants(barcode,
                              ham_dist=0):
    """Generate all barcode variants within a given Hamming distance.
    
    Args:
        barcode (str): original barcode sequence (e.g., 'ACGT').
        ham_dist (int): maximum Hamming distance (default 0; adjusted to 1 if 0).
    
    Returns:
        list: List of unique barcode variants, including the original.
    """
    # Use 'N' only if ham_dist is 0; otherwise, full nucleotide set
    nucleotides = ['N'] if ham_dist == 0 else ['A', 'C', 'G', 'T', 'N']
    max_dist = max(1, ham_dist)  # Ensure at least distance 1
    
    # makeing a set of barocde and keep original barcode
    variants = {barcode}
    
    # iterate over each Hamming distance up to max_dist
    for dist in range(1, max_dist + 1):
        # Choose 'dist' positions to mutate
        for pos_combo in combinations(range(len(barcode)), dist):
            # Generate all possible nucleotide combinations for these positions
            for sub_combo in product(nucleotides, repeat=dist):
                # Create variant by substituting at chosen positions
                new_barcode = list(barcode)
                for pos, nu in zip(pos_combo, sub_combo):
                    new_barcode[pos] = nu
                variants.add(''.join(new_barcode))
    
    return list(variants)


# --------------------------------------------------
def load_barcode_df(file_path,
                    ham_dist=0,
                    sep=None,
                    names=['barcode', 'cellID'], 
                    index_col=0,
                    compatible_bc='compatible_bcs',
                    ):
    """Load a barcode DataFrame from a file and compute compatible barcodes.

    Args:
        file_path (str): path to the tsv: contain barcode and cellID/condition
        ham_dist: Hamming distance
        sep (str, optional): separator (',' or '\t'). If None, tries to detect automatically
        names (list): column names setting. Defaults to ['bc', 'cellID']
        index_col (int): index column. Defaults to 0 (first column)

    Returns:
        pandas.DataFrame: DataFrame with barcode data and compatible barcodes.
    """
    
    # If sep is not provided, try to detect it
    if sep is None:
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
            if '\t' in first_line:
                sep = '\t'
            else:
                sep = ','  # Default to comma if tab not found
    
    # Read the file with the specified or detected separator
    barcode_df = pd.read_csv(file_path, sep=sep, names=names, index_col=index_col)
    
    # Add compatible barcodes column using the provided Hamming distance from args
    barcode_df[compatible_bc] = barcode_df.apply(lambda x: generate_barcode_variants(x.name, ham_dist), axis=1)
    
    return barcode_df


# --------------------------------------------------
def extract_unique_variants(barcode_df, 
                            condition='cellID',
                            compatible_bc='compatible_bcs'):
    """
    filter compatible barcode variants that appear once for each cellID
    
    Args:
        bc_df (DataFrame): 'compatible_bcs' and condition as 'cellID' 
    
    Returns:
        DataFrame with unique barcodes as index and their 'cellID' and 'original' index as columns
    """
    # count uniqueness of compatible_bcs across barcodes
    count_all_comp_bc = Counter([comp_bc for barcode in barcode_df.index for comp_bc in barcode_df.loc[barcode, compatible_bc]])
    
    # loop the compatible barcode:if appear once across set of compatible barcode pool,
    # then assign with original barcode
    unique_bc = {
        comp_bc: {condition: row[condition], 'original': barcode}
        for barcode, row in barcode_df.iterrows()
        for comp_bc in row[compatible_bc]
        if count_all_comp_bc[comp_bc] == 1
    }
    
    # convert to dataframe
    unique_bc_df = pd.DataFrame(unique_bc).T
    
    return unique_bc_df

