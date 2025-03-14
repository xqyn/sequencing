
print('Importing barcode_utils.py')

from itertools import product, combinations

def generate_barcode_variants(barcode, ham_dist=0):
    """ 
    Returns list of barcode sequences within given Hamming distance.
    """
    # set possible nucleotides and adjust ham_dist
    nucleotides = ['N'] if ham_dist == 0 else ['A', 'C', 'G', 'T', 'N']
    ham_dist = 1 if ham_dist == 0 else ham_dist
    
    # keep original barcode in barcode set
    barcode_variant = {barcode}
    
    # loop over all possible hamming distances
    for distance in range(1, ham_dist + 1):
        # generate all possible positions     
        positions = combinations(range(len(barcode)), distance)
        
        # for each position set
        for pos_set in positions:
            # try each nucleotide substitution
            subs = product(nucleotides, repeat=distance)
            # apply substitutions
            for sub in subs:
                new_barcode = list(barcode)
                for p, n in zip(pos_set, sub):
                    new_barcode[p] = n
                barcode_variant.add(''.join(new_barcode))
    return list(barcode_variant)
