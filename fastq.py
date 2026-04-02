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
from dataclasses import dataclass


# --- configure logging --------------------------------------------------
import logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)


# --- helper --------------------------------------------------
def ascii_to_phred(quality_string, offset=33):
    """Convert an ASCII quality to Phred scores"""
    return [ord(char) - offset for char in quality_string]


def phred_to_ascii(phred_scores: list) -> str:
    """Convert integer PHRED scores to ASCII using PHRED+33 encoding."""
    return ''.join(chr(q + 33) for q in phred_scores)


def complement(seq, reverse=False):
    """Return the complement of a sequence"""
    mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', '-':'-',
               'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}
    complemented = ''.join([mapping.get(base, base) for base in seq])
    return complemented[::-1] if reverse else complemented


# --- class --------------------------------------------------
@dataclass(frozen=True)
class Fastq:
    """
    Define a fastq class for storing FASTQ records.
    
    Attributes:
        iden    (str): Identifier of the sequence (e.g., read name).
        seq     (str): DNA sequence string.
        qual    (str): Quality string (ASCII-encoded scores).
        phred   (Optional[int]): List of Phred quality scores derived from qual.
        length  (int): Length of the sequence.
    """
    iden:   str     # identifier
    seq:    str     # sequence
    qual:   str     # quality
    phred:  Optional[List[int]]
    length: int     # length of the sequence
        

# --- SeqQuery --------------------------------------------------
def SeqQuery(fastq_filename: str,
             phred: bool = True,
             lines: int = None) -> Iterator[Fastq]:
    """
    Parse a FASTQ (or gzipped FASTQ) file and yield Fastq records.
    Args:
        fastq_filename  (str)   : path to the FASTQ file (.fastq or .fastq.gz).
        phred           (bool)  : whether to compute Phred scores (default: True).
        lines           (int)   : number of FASTQ records to yield (default: None = all).
    
    Yields:
        fastq: A fastq object for each record.
    
    Raises:
        ValueError: If FASTQ records are incomplete or sequence/quality lengths differ.
    """
    opener = gzip.open if fastq_filename.endswith('.gz') else open
    
    with opener(fastq_filename, 'rt') as fh:
        records = zip(*[iter(fh)] * 4)
        
        for i, (header, seq, plus, qual) in enumerate(islice(records, lines)):
            header = header.strip()
            seq    = seq.strip()
            plus   = plus.strip()
            qual   = qual.strip()
            
            # 1. Header must start with '@'
            if not header.startswith('@'):
                raise ValueError(f"Record {i}: header must start with '@', got: {header!r}")
            
            # 2. Sequence must not be empty
            if not seq:
                raise ValueError(f"Record {i}: empty sequence for read '{header}'")
            
            # 3. '+' separator line must start with '+'
            if not plus.startswith('+'):
                raise ValueError(f"Record {i}: expected '+' separator, got: {plus!r}")
            
            # 4. Quality must not be empty
            if not qual:
                raise ValueError(f"Record {i}: empty quality string for read '{header}'")
            
            # 5. Sequence and quality lengths must match
            if len(seq) != len(qual):
                raise ValueError(
                    f"Record {i}: sequence length ({len(seq)}) != "
                    f"quality length ({len(qual)}) for read '{header}'"
                )
            
            # 6. Quality scores must be in valid Phred+33 ASCII range (33–126)
            if not all(33 <= ord(c) <= 126 for c in qual):
                raise ValueError(f"Record {i}: quality string contains out-of-range ASCII for read '{header}'")
            
            yield Fastq(
                iden=header,
                seq=seq,
                qual=qual,
                phred=ascii_to_phred(qual) if phred else None,
                length=len(seq)
            )
        
        # 7. check for trailing incomplete record (truncated file)
        leftover = list(islice(fh, 3))
        if leftover:
            raise ValueError(
                f"File ended with an incomplete FASTQ record "
                f"({len(leftover)} trailing line(s))"
            )



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