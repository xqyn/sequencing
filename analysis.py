#!/bin/bash python

'''
BELT - XQ Ng - 14 March 2025
Leiden UMC
'''

# import --------------------------------------------------
import sys
custom_paths = ['./utils/']

for p in custom_paths:
    if p not in sys.path:
        sys.path.append(p)

from fastq import *
# --------------------------------------------------

fastq_R1_name='data/ABE_0d0_R1.fastq.gz'
# iteretation
fastq_R1 = fastq_seq(fastq_R1_name)
# load into a list: cause memory
fastq_R1_list=list(fastq_seq(fastq_R1_name))


fastq_R2_name='data/ABE_0d0_R2.fastq.gz'
fastq_R2=list(fastq_seq(fastq_R2_name))

fastq_filename='abe_R1.txt'

test=fastq_seq(fastq_filename)


# --------------------------------------------------
# Assuming fastq_R1[1] is your fastq object
fastq_record = next(fastq_R1)

# Extract each element
identifier = fastq_record.iden  # Gets the identifier string
sequence = fastq_record.seq    # Gets the sequence string
quality = fastq_record.qual    # Gets the quality string
phred_scores = fastq_record.phred  # Gets the list of Phred scores
seq_length = fastq_record.length   # Gets the sequence length

# Print or use the extracted elements
print("Identifier:", identifier)
print("Sequence:", sequence)
print("Quality:", quality)
print("Phred Scores:", phred_scores)
print("Length:", seq_length)


sequences = [record.seq for record in fastq_R1]