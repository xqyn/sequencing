#!/bin/bash python
import subprocess as sb
import gzip

fastq_filename='data/ABE_0d0_R1.fastq.gz'

# --------------------------------------------------
# class
class fastq(NamedTuple):
    
    """
    Define a fastq class for storing fastq records: 
    iden:   identifier
    seq:    sequence
    qual:   quality
    """
    iden: str   # identifier
    seq: str    # sequence
    qual: str   # quality


# --------------------------------------------------
def fastq_seq(fastq_filename):
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
            if seq_num % 4 == 3: # for evey 4th line
                seq_record = fastq(
                    iden=temp_seq[0],
                    seq=temp_seq[1],
                    qual=temp_seq[3]
                )
                fastq_list.append(seq_record)
    return fastq_list



def seq_len(fastq_list):
    """
    Giving a fastq_list containing fastq.seq
    return the lengths of the sequences
    """
    
    return [len(sequencing.seq) for sequencing in fastq_list]