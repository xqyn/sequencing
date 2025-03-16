
print('Importing barcode_utils.py')

from itertools import combinations, product

def generate_barcode_variants(barcode, ham_dist=0):
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
