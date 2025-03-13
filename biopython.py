from Bio import SeqIO

# Path to the FASTQ file
fastq_file = "abe_R1.txt"

# Parse the FASTQ file and extract quality scores
for record in SeqIO.parse(fastq_file, "fastq"):
    print(f"Sequence ID: {record.id}")
    print(f"Sequence: {record.seq}")
    # Quality scores are automatically converted to Phred values (numbers)
    quality_scores = record.letter_annotations["phred_quality"]
    print(f"Quality Scores: {quality_scores}")