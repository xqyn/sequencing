from Bio import SeqIO
from Bio.Seq import Seq

# --------------------------------------------------
# turorial
for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

# --------------------------------------------------
# Sequence objects

## Sequences act like strings
my_seq = Seq("GATCG")
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))

print(my_seq[0])
"AAAA".count("AA")
Seq("AAAA").count("AA")
my_seq.count('G')

my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
my_seq.count("G")

# count percentage of G and C in the sequence
round((my_seq.count('G') + my_seq.count('C'))*100 / len(my_seq),2)

from Bio.SeqUtils import gc_fraction
gc_fraction(my_seq)

## Slicing a sequence
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
my_seq[4:12]


## reverse string
my_seq[::-1]

## Turning Seq objects into strings
str(my_seq)

## Concatenating or adding sequences

list_of_seqs = [Seq("ACGT"), Seq("AACC"), Seq("GGTT")]
concatenated = Seq("")
for s in list_of_seqs:
    concatenated += s

### or use join
contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTGCA")]
spacer = Seq("N" * 10)
spacer.join(contigs)
Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTGCA')


## Changing case
dna_seq = Seq("acgtACGT")
dna_seq.upper()
dna_seq.lower()
"GTAC" in dna_seq.upper()

# --------------------------------------------------
# Sequence annotation objects


# --------------------------------------------------
# Sequence Input/Output
from Bio import SeqIO
help(SeqIO)

## Parsing or Reading Sequences

from Bio import SeqIO

for record in SeqIO.parse("m_cold.fasta", "fasta"):
    print(record.id, len(record))