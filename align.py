def dna_alignment(seq1, seq2):
    """
    Function to align two DNA sequences and calculate similarity
    Parameters:
        seq1: First DNA sequence (string)
        seq2: Second DNA sequence (string)
    Returns:
        alignment visualization and similarity score
    """
    # Ensure sequences are uppercase
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    
    # Get lengths of sequences
    len1, len2 = len(seq1), len(seq2)
    
    # If sequences are of different lengths, pad the shorter one with gaps
    if len1 < len2:
        seq1 = seq1 + '-' * (len2 - len1)
    elif len2 < len1:
        seq2 = seq2 + '-' * (len1 - len2)
    
    # Scoring system
    match_score = 1
    mismatch_score = -1
    gap_penalty = -2
    
    # Calculate similarity
    score = 0
    matches = 0
    mismatches = 0
    gaps = 0
    
    # Create alignment visualization
    alignment = ""
    
    for i in range(max(len1, len2)):
        if seq1[i] == seq2[i]:
            if seq1[i] != '-':
                score += match_score
                matches += 1
                alignment += "|"
            else:
                score += gap_penalty
                gaps += 1
                alignment += " "
        else:
            score += mismatch_score
            mismatches += 1
            alignment += "."
    
    # Calculate similarity percentage
    total_positions = max(len1, len2)
    similarity_percent = (matches / total_positions) * 100 if total_positions > 0 else 0
    
    # Print results
    print(f"Sequence 1: {seq1}")
    print(f"Alignment : {alignment}")
    print(f"Sequence 2: {seq2}")
    print(f"\nAlignment Statistics:")
    print(f"Total length: {total_positions}")
    print(f"Matches: {matches}")
    print(f"Mismatches: {mismatches}")
    print(f"Gaps: {gaps}")
    print(f"Score: {score}")
    print(f"Similarity: {similarity_percent:.2f}%")
    
    return score, similarity_percent

# Example usage
def main():
    # Test sequences
    sequence1 = "ATCGATCG"
    sequence2 = "ATGGATGG"
    
    print("DNA Sequence Alignment Results:")
    print("-" * 30)
    
    score, similarity = dna_alignment(sequence1, sequence2)

if __name__ == "__main__":
    main()