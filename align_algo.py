# --------------------------------------------------
# needleman_wunsch

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    # Initialize the scoring matrix with zeros
    rows, cols = len(seq1) + 1, len(seq2) + 1
    matrix = [[0 for _ in range(cols)] for _ in range(rows)]
    
    # Fill first row and column with gap penalties
    for i in range(rows):
        matrix[i][0] = i * gap
    for j in range(cols):
        matrix[0][j] = j * gap
    
    # Fill the scoring matrix
    for i in range(1, rows):
        for j in range(1, cols):
            match_score = match if seq1[i-1] == seq2[j-1] else mismatch
            diagonal = matrix[i-1][j-1] + match_score  # Match/mismatch
            vertical = matrix[i-1][j] + gap            # Gap in seq2
            horizontal = matrix[i][j-1] + gap          # Gap in seq1
            matrix[i][j] = max(diagonal, vertical, horizontal)
    
    # Traceback to get the aligned sequences
    aligned1, aligned2 = [], []
    i, j = len(seq1), len(seq2)
    while i > 0 or j > 0:
        if i > 0 and j > 0 and matrix[i][j] == matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif i > 0 and matrix[i][j] == matrix[i-1][j] + gap:
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            i -= 1
        else:
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            j -= 1
    
    # Reverse the sequences since we built them backwards
    aligned1 = ''.join(reversed(aligned1))
    aligned2 = ''.join(reversed(aligned2))
    
    return aligned1, aligned2, matrix[len(seq1)][len(seq2)]

# Example usage
seq1 = "CAT"
seq2 = "CGT"
aligned1, aligned2, score = needleman_wunsch(seq1, seq2)
print(f"Aligned seq1: {aligned1}")
print(f"Aligned seq2: {aligned2}")
print(f"Score: {score}")

# --------------------------------------------------
def smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-2):
    # Initialize the scoring matrix with zeros
    rows, cols = len(seq1) + 1, len(seq2) + 1
    matrix = [[0 for _ in range(cols)] for _ in range(rows)]
    
    # Track the maximum score and its position for traceback
    max_score = 0
    max_pos = (0, 0)
    
    # Fill the scoring matrix
    for i in range(1, rows):
        for j in range(1, cols):
            match_score = match if seq1[i-1] == seq2[j-1] else mismatch
            diagonal = matrix[i-1][j-1] + match_score  # Match/mismatch
            vertical = matrix[i-1][j] + gap            # Gap in seq2
            horizontal = matrix[i][j-1] + gap          # Gap in seq1
            matrix[i][j] = max(0, diagonal, vertical, horizontal)  # Key difference: reset to 0 if negative
            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_pos = (i, j)
    
    # Traceback from the max score position
    aligned1, aligned2 = [], []
    i, j = max_pos
    while i > 0 and j > 0 and matrix[i][j] > 0:
        if matrix[i][j] == matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif matrix[i][j] == matrix[i-1][j] + gap:
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            i -= 1
        else:
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            j -= 1
    
    # Reverse the sequences
    aligned1 = ''.join(reversed(aligned1))
    aligned2 = ''.join(reversed(aligned2))
    
    return aligned1, aligned2, max_score

# Example usage
seq1 = "TGTTACGG"
seq2 = "GGTTGACTA"
aligned1, aligned2, score = smith_waterman(seq1, seq2)
print(f"Aligned seq1: {aligned1}")
print(f"Aligned seq2: {aligned2}")
print(f"Score: {score}")



# --------------------------------------------------
# combine both method

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    # Same as before, but we’ll return the matrix too for region extraction
    rows, cols = len(seq1) + 1, len(seq2) + 1
    matrix = [[0 for _ in range(cols)] for _ in range(rows)]
    for i in range(rows):
        matrix[i][0] = i * gap
    for j in range(cols):
        matrix[0][j] = j * gap
    for i in range(1, rows):
        for j in range(1, cols):
            match_score = match if seq1[i-1] == seq2[j-1] else mismatch
            diagonal = matrix[i-1][j-1] + match_score
            vertical = matrix[i-1][j] + gap
            horizontal = matrix[i][j-1] + gap
            matrix[i][j] = max(diagonal, vertical, horizontal)
    
    # Traceback
    aligned1, aligned2 = [], []
    i, j = len(seq1), len(seq2)
    while i > 0 or j > 0:
        if i > 0 and j > 0 and matrix[i][j] == matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif i > 0 and matrix[i][j] == matrix[i-1][j] + gap:
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            i -= 1
        else:
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            j -= 1
    aligned1 = ''.join(reversed(aligned1))
    aligned2 = ''.join(reversed(aligned2))
    return aligned1, aligned2, matrix[len(seq1)][len(seq2)]

def smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-2):
    rows, cols = len(seq1) + 1, len(seq2) + 1
    matrix = [[0 for _ in range(cols)] for _ in range(rows)]
    max_score, max_pos = 0, (0, 0)
    for i in range(1, rows):
        for j in range(1, cols):
            match_score = match if seq1[i-1] == seq2[j-1] else mismatch
            diagonal = matrix[i-1][j-1] + match_score
            vertical = matrix[i-1][j] + gap
            horizontal = matrix[i][j-1] + gap
            matrix[i][j] = max(0, diagonal, vertical, horizontal)
            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_pos = (i, j)
    
    # Traceback
    aligned1, aligned2 = [], []
    i, j = max_pos
    while i > 0 and j > 0 and matrix[i][j] > 0:
        if matrix[i][j] == matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif matrix[i][j] == matrix[i-1][j] + gap:
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            i -= 1
        else:
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            j -= 1
    aligned1 = ''.join(reversed(aligned1))
    aligned2 = ''.join(reversed(aligned2))
    return aligned1, aligned2, max_score

# Hybrid approach
def hybrid_alignment(short_seq, long_seq):
    # Step 1: Global alignment with Needleman-Wunsch
    aligned_short, aligned_long, global_score = needleman_wunsch(short_seq, long_seq)
    print("Global Alignment:")
    print(f"Short: {aligned_short}")
    print(f"Long:  {aligned_long}")
    print(f"Global Score: {global_score}")
    
    # Extract the region of the long sequence aligned to the short sequence
    # Remove gaps from aligned_short to find its span in aligned_long
    short_no_gaps = ''.join([c for c in aligned_short if c != '-'])
    start_idx = aligned_long.index(short_no_gaps[0])  # First matching char
    end_idx = start_idx + len(short_no_gaps)  # Span of short seq
    long_region = long_seq[start_idx:end_idx]  # Extract the region
    
    print(f"\nExtracted region from long sequence: {long_region}")
    
    # Step 2: Local alignment with Smith-Waterman on the extracted region
    local_aligned_short, local_aligned_region, local_score = smith_waterman(short_seq, long_region)
    print("\nLocal Alignment (Refined):")
    print(f"Short: {local_aligned_short}")
    print(f"Region: {local_aligned_region}")
    print(f"Local Score: {local_score}")

# Test it
short_seq = "GTTAC"
long_seq = "TGATTGTTACGG"
hybrid_alignment(short_seq, long_seq)