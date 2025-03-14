def pairwise_align(seq1, seq2, mode="global", match_score=1, mismatch_score=-1, open_gap_score=-2, extend_gap_score=-1):
    """
    Mimics Biopython's Bio.Align.PairwiseAligner for global/local alignment with affine gaps.
    Args:
        seq1, seq2: Strings to align (e.g., "TACCG", "ACG").
        mode: "global" (Needleman-Wunsch) or "local" (Smith-Waterman).
        match_score, mismatch_score: Scores for match/mismatch.
        open_gap_score, extend_gap_score: Affine gap penalties (negative values).
    Returns:
        aligned_seq1, aligned_seq2 (strings), score (float).
    """
    # Initialize matrices: score matrix (F) and gap matrices (Ix, Iy) for affine gaps
    rows, cols = len(seq1) + 1, len(seq2) + 1
    F = [[0] * cols for _ in range(rows)]  # Main scoring matrix
    Ix = [[float('-inf')] * cols for _ in range(rows)]  # Gaps in seq1 (horizontal)
    Iy = [[float('-inf')] * cols for _ in range(rows)]  # Gaps in seq2 (vertical)

    # For local alignment, track the maximum score and its position
    max_score = 0
    max_pos = (0, 0)

    # Initialize first row and column (global mode only)
    if mode == "global":
        for i in range(rows):
            F[i][0] = open_gap_score + (i - 1) * extend_gap_score if i > 0 else 0
            Iy[i][0] = F[i][0]  # Gap in seq2
        for j in range(cols):
            F[0][j] = open_gap_score + (j - 1) * extend_gap_score if j > 0 else 0
            Ix[0][j] = F[0][j]  # Gap in seq1

    # Fill the matrices using affine gap logic
    for i in range(1, rows):
        for j in range(1, cols):
            # Match or mismatch score
            match = match_score if seq1[i-1] == seq2[j-1] else mismatch_score

            # Affine gap options
            Ix[i][j] = max(
                F[i][j-1] + open_gap_score,    # Open gap in seq1
                Ix[i][j-1] + extend_gap_score  # Extend gap in seq1
            )
            Iy[i][j] = max(
                F[i-1][j] + open_gap_score,    # Open gap in seq2
                Iy[i-1][j] + extend_gap_score  # Extend gap in seq2
            )
            F[i][j] = max(
                F[i-1][j-1] + match,  # Diagonal: match/mismatch
                Ix[i][j],             # Gap in seq1
                Iy[i][j]              # Gap in seq2
            )

            # Local alignment: reset to 0 if negative, track max score
            if mode == "local":
                F[i][j] = max(0, F[i][j])
                if F[i][j] > max_score:
                    max_score = F[i][j]
                    max_pos = (i, j)

    # Traceback
    aligned1, aligned2 = [], []
    if mode == "global":
        i, j = len(seq1), len(seq2)  # Start from bottom-right
        score = F[i][j]
    else:  # Local
        i, j = max_pos  # Start from max score position
        score = max_score

    while (mode == "global" and (i > 0 or j > 0)) or (mode == "local" and F[i][j] > 0):
        if i > 0 and j > 0 and F[i][j] == F[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score):
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif i > 0 and F[i][j] == Iy[i][j]:  # Gap in seq2
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            i -= 1
        elif j > 0 and F[i][j] == Ix[i][j]:  # Gap in seq1
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            j -= 1
        else:  # Handle edge cases or ambiguities
            break

    # Reverse the aligned sequences
    aligned1 = ''.join(reversed(aligned1))
    aligned2 = ''.join(reversed(aligned2))

    return aligned1, aligned2, score

# Example usage
seq1 = "TACCG"
seq2 = "ACG"

# Global alignment (Needleman-Wunsch with affine gaps)
aligned1_g, aligned2_g, score_g = pairwise_align(seq1, seq2, mode="global")
print("Global Alignment:")
print(f"{aligned1_g}\n{aligned2_g}\nScore: {score_g}")

# Local alignment (Smith-Waterman with affine gaps)
aligned1_l, aligned2_l, score_l = pairwise_align(seq1, seq2, mode="local")
print("\nLocal Alignment:")
print(f"{aligned1_l}\n{aligned2_l}\nScore: {score_l}")