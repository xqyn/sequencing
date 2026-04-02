'''
project: belt
XQ - Leiden UMC
FASTQ modifier
update:
    - juli 23: add extract_pair_fastq, 
    - juli 24-25: add FastqPair 
    - aug 6: write_to_fastq
    - aug 9: add write_fastq_queries
'''

import gzip
from typing import List, Tuple, Optional, NamedTuple, Iterator
from itertools import islice
from Bio import SeqIO
from sequence import ascii_to_phred, phred_to_ascii, complement

# configure logging
import logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)


# function ---
def seq_lenght(fastq_list):
    """Calculate the lengths of sequences"""
    return [len(sequencing.seq) for sequencing in fastq_list]

# def ascii_to_phred(quality_string, offset=33):
#     """Convert an ASCII quality to Phred scores"""
#     return [ord(char) - offset for char in quality_string]


# class --------------------------------------------------
class Fastq(NamedTuple):
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
        

# SeqQuery --------------------------------------------------
def SeqQuery(fastq_filename: str,
             phred_offset: int = 33,
             phred: bool = True,
             lines: int = None) -> Iterator[Fastq]:
    """
    Args:
        fastq_filename: The name of the fastq file.
        phred_offset: Phred score offset for quality conversion (default: 33).
        lines: Number of FASTQ records to yield (default: None, yields all records).
    
    Yields:
        fastq: A fastq object for each record.
    
    Raises:
        ValueError: If FASTQ records are incomplete or sequence/quality lengths differ.
    """
    
    fastq_list = []
    opener = gzip.open if fastq_filename.endswith('.gz') else open
    with opener(fastq_filename, 'rt') as sequencing:
        # Set iterator with optional limit using islice
        seq_iterator = islice(sequencing, None) if lines is None else islice(sequencing, lines * 4)
        
        # seting objects for temp sequecing
        temp_seq = [None, None, None, None]
        for seq_num, line in enumerate(seq_iterator):
            temp_seq[seq_num % 4] = line.strip()    # store line in temp_seq
            if seq_num % 4 == 3:                    # for evey 4th line
                if None in temp_seq:
                    raise ValueError(f"Incomplete FASTQ record at line {seq_num - 2}")
                if len(temp_seq[1]) != len(temp_seq[3]):
                    raise ValueError(f"Sequence and quality lengths differ at record starting line {seq_num - 2}")
                seq_record = Fastq(
                    iden=temp_seq[0],
                    seq=temp_seq[1],
                    qual=temp_seq[3],
                    phred=ascii_to_phred(temp_seq[3]) if phred else None,
                    length=len(temp_seq[3])
                )
                yield seq_record
    # Check for trailing incomplete records
        if seq_num % 4 != 3:
            remaining_lines = (seq_num + 1) % 4
            if remaining_lines != 0:
                raise ValueError(f"File ended with incomplete FASTQ record (lines: {seq_num + 1})")



# ---Extract pair
class FastqPair(NamedTuple):
    """Container for paired-end FASTQ record data."""
    read_id: str            # Read ID DNA 
    r1_seq: str             # Read 1 DNA sequence
    r1_phred: List[int]     # Read 1 Phred quality scores
    r2_seq: str             # Read 2 DNA sequence
    r2_phred: List[int]     # Read 2 Phred quality scores
    r2_seq_rev: str         # Reverse complement of Read 2 sequence
    r2_phred_rev: List[int] # Reversed Read 2 Phred scores
    


# --- based on query_id
def read_id_fastq(read_id,
                  filepath: str) -> Tuple[Optional[str], Optional[List[int]]]:
    """Read a single FASTQ file and extract sequence and quality."""
    opener = gzip.open if filepath.endswith('.gz') else open
    
    try:
        with opener(filepath, 'rt') as file:
            for record in SeqIO.parse(file, 'fastq'):
                if record.id == read_id[1:]:
                    return str(record.seq), record.letter_annotations["phred_quality"]
    except FileNotFoundError:
        logger.error(f"Error: File {filepath}: not found")
    except Exception as e:
        logger.error(f"Error reading {filepath}: {str(e)}")
        
    return None, None

def read_id_SeqQuery(read_id,
                  filepath: str) -> Tuple[Optional[str], Optional[List[int]]]:
    """Read a single FASTQ file and extract sequence and quality."""
    fastq = SeqQuery(filepath, phred=False)
    try:
        for query in fastq:
            if query.iden == read_id:
                return str(query.seq), ascii_to_phred(query.qual)
    except FileNotFoundError:
        logger.error(f"Error: File {filepath}: not found")
    except Exception as e:
        logger.error(f"Error reading {filepath}: {str(e)}")
        
    return None, None


def extract_pair_fastq(read_id: str, 
                       read1_file: str, 
                       read2_file: str) -> Optional[FastqPair]:
    """
    Extract paired-end FASTQ records for a given read ID.
    
    Args:
        read_id (str): Identifier of the sequence (with or without '@').
        read1_file (str): Path to the Read 1 FASTQ file (gzipped or plain).
        read2_file (str): Path to the Read 2 FASTQ file (gzipped or plain).
    
    Returns:
        Record or None if read not found or error occurs.
    """
    #logger.info(f"Processing FASTQ ID: {read_id}")
    # Ensure read_id starts with '@' for FASTQ format
    read_id = read_id if read_id.startswith('@') else '@' + read_id
    
    r1_seq, r1_phred, r2_seq, r2_phred = None, None, None, None
    
    if read1_file:
        #r1_seq, r1_phred = read_id_fastq(read_id, read1_file)
        r1_seq, r1_phred = read_id_SeqQuery(read_id, read1_file)
    if read2_file:
        #r2_seq, r2_phred = read_id_fastq(read_id, read2_file)
        r2_seq, r2_phred = read_id_SeqQuery(read_id, read2_file)
    
    if not r1_seq and read1_file:
        print(f"Error: Read {read_id} not found in Read 1 file.")
        return None
    if not r2_seq and read2_file:
        print(f"Error: Read {read_id} not found in Read 2 file.")
        return None
    
    # Generate reverse complement and reversed quality scores
    r2_seq_rev = complement(r2_seq, reverse=True) if r2_seq else None
    r2_phred_rev = r2_phred[::-1] if r2_phred else None
    
    return FastqPair(
        read_id=read_id,
        r1_seq=r1_seq,
        r1_phred=r1_phred,
        r2_seq=r2_seq,
        r2_phred=r2_phred,
        r2_seq_rev=r2_seq_rev,
        r2_phred_rev=r2_phred_rev
    )
    
# # write out --------------------------------------------------
def write_to_fastq(output_path, query_id, sequence, seq_score, phred_str=True):
    """
    Write sequence pair to FASTQ.gz file.
    """
    with gzip.open(output_path, 'at', encoding='utf-8') as f:
        quality = phred_to_ascii(seq_score) if phred_str else seq_score
        f.write(f"{query_id}\n{sequence}\n+\n{quality}\n")
        
        
# def write_to_fastq(output_path, query_id, sequence, seq_score, phred_str=True):
#     """
#     Write sequence pair to FASTQ.gz file.
#     """
#     # Ensure inputs are strings
#     query_id = query_id.decode('utf-8') if isinstance(query_id, bytes) else query_id
#     sequence = sequence.decode('utf-8') if isinstance(sequence, bytes) else sequence
#     seq_score = seq_score.decode('utf-8') if isinstance(seq_score, bytes) else seq_score
#     quality = phred_to_ascii(seq_score) if phred_str else seq_score
#     if isinstance(quality, bytes):
#         quality = quality.decode('utf-8')
#         logger.warning(f"phred_to_ascii returned bytes for query_id {query_id}")
    
#     # Encode FASTQ text to bytes and write in binary mode
#     fastq_entry = f"@{query_id}\n{sequence}\n+\n{quality}\n".encode('utf-8')
#     with gzip.open(output_path, 'ab') as f:
#         f.write(fastq_entry)

def write_fastq_queries(query_list, temp_file, phred_str=True):
    with gzip.open(temp_file, 'wt') as f:
        for query_id, sequence, seq_score in query_list:
            # identifier
            f.write(f'@{query_id}\n')
            # sequence
            f.write(f'{sequence}\n')
            # line
            f.write('+\n')
            # convert Phred scores to ASCII (Phred+33) and write
            quality = phred_to_ascii(seq_score) if phred_str else seq_score
            f.write(f'{quality}\n')