#!/bin/bash python

'''
BELT - XQ Ng - 14 March 2025
Leiden UMC
'''

import gzip
import subprocess as sb
from typing import NamedTuple, List, Iterator


# --------------------------------------------------
# class
class fastq(NamedTuple):
    """
    Define a fastq class for storing FASTQ records.

    Attributes:
        iden (str): Identifier of the sequence (e.g., read name).
        seq (str): DNA sequence string.
        qual (str): Quality string (ASCII-encoded scores).
        phred (List[int]): List of Phred quality scores derived from qual.
        length (int): Length of the sequence.
    """
    iden:   str     # identifier
    seq:    str     # sequence
    qual:   str     # quality
    phred:  List[int]     # phred quality score 
    length: int     # length of the sequence
        
    


# --------------------------------------------------
def fastq_seq(fastq_filename: str, 
              phred_offset: int = 33) -> Iterator[fastq]:
    """
    Read a fastq file and return a list of fastq objects.
    fastq_filename: the name of the fastq file
    """
    
    fastq_list = []
    opener = gzip.open if fastq_filename.endswith('.gz') else open
    with opener(fastq_filename, 'rt') as sequencing:
        # seting objects for temp sequecing
        temp_seq = [None, None, None, None]
        for seq_num, line in enumerate(sequencing):
            temp_seq[seq_num % 4] = line.strip()    # store line in temp_seq
            if seq_num % 4 == 3:                    # for evey 4th line
                if None in temp_seq:
                    raise ValueError(f"Incomplete FASTQ record at line {seq_num - 2}")
                if len(temp_seq[1]) != len(temp_seq[3]):
                    raise ValueError(f"Sequence and quality lengths differ at record starting line {seq_num - 2}")
                seq_record = fastq(
                    iden=temp_seq[0],
                    seq=temp_seq[1],
                    qual=temp_seq[3],
                    phred=ascii_to_phred(temp_seq[3]),
                    length=len(temp_seq[3])
                )
                yield seq_record
    # Check for trailing incomplete records
        if seq_num % 4 != 3:
            remaining_lines = (seq_num + 1) % 4
            if remaining_lines != 0:
                raise ValueError(f"File ended with incomplete FASTQ record (lines: {seq_num + 1})")

def seq_lenght(fastq_list):
    """Calculate the lengths of sequences"""
    return [len(sequencing.seq) for sequencing in fastq_list]

def ascii_to_phred(quality_string, offset=33):
    """Convert an ASCII quality to Phred scores"""
    return [ord(char) - offset for char in quality_string]