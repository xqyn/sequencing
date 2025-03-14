def find_indels_substitutions(read_seq_al, ref_seq_al, _include_indx):
    ref_positions = []
    substitutions = []
    deletions = []
    insertions = []
    
    include_indx_set = set(_include_indx)
    idx = 0
    start_deletion = -1
    start_insertion = -1
    insertion_size = 0
    
    for i, (ref_c, read_c) in enumerate(zip(ref_seq_al, read_seq_al)):
        if ref_c != '-':
            ref_positions.append(idx)
            if read_c != ref_c and read_c != '-' and read_c != 'N':
                substitutions.append((idx, read_c))
            if start_insertion != -1:
                insertions.append((start_insertion, idx, insertion_size))
                start_insertion = -1
                insertion_size = 0
            idx += 1
        else:
            ref_positions.append(-idx if idx > 0 else -1)
            if idx > 0 and start_insertion == -1:
                start_insertion = idx - 1
            insertion_size += 1
        
        if read_c == '-' and start_deletion == -1:
            start_deletion = ref_positions[i] if i > 0 else 0
        elif read_c != '-' and start_deletion != -1:
            deletions.append((start_deletion, ref_positions[i]))
            start_deletion = -1
    
    if start_deletion != -1:
        deletions.append((start_deletion, ref_positions[-1]))
    
    filtered_substitutions = [(pos, val) for pos, val in substitutions if pos in include_indx_set]
    filtered_deletions = [(start, end) for start, end in deletions 
                         if include_indx_set.intersection(range(start, end))]
    filtered_insertions = [(start, end, size) for start, end, size in insertions 
                          if start in include_indx_set and end in include_indx_set]
    
    substitution_n = len(filtered_substitutions)
    deletion_n = sum(end - start for start, end in filtered_deletions)
    insertion_n = sum(size for _, _, size in filtered_insertions)
    
    return {
        "all_substitution_positions": [pos for pos, _ in substitutions],
        "substitution_positions": [pos for pos, _ in filtered_substitutions],
        "all_substitution_values": [val for _, val in substitutions],
        "substitution_values": [val for _, val in filtered_substitutions],
        "substitution_n": substitution_n,
        "all_deletion_coordinates": deletions,
        "deletion_coordinates": filtered_deletions,
        "deletion_sizes": [end - start for start, end in filtered_deletions],
        "deletion_n": deletion_n,
        "all_insertion_coordinates": [(start, end) for start, end, _ in insertions],
        "insertion_coordinates": [(start, end) for start, end, _ in filtered_insertions],
        "insertion_sizes": [size for _, _, size in filtered_insertions],
        "insertion_n": insertion_n,
        "ref_positions": ref_positions
    }
    
ref_seq_al = "AATT-GCC"
read_seq_al = "A-TTTGGC"
_include_indx = [1, 2, 3]

result = find_indels_substitutions(read_seq_al, ref_seq_al, _include_indx)
print(result)