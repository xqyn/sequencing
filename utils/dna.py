
#!/bin/bash python

'''
BELT - XQ Ng - 14 March 2025
Leiden UMC
''' 
print('Importing dna_utils.py')

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