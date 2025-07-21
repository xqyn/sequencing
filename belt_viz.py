"""
project: belt
XQ - Leiden UMC
visualization of squencing alignments
update: 2025
    - juni 16: plot_base_grid
    - juni 26: fix x_label for position base + fsize adds
    - juli 15: cigar_recover (checking docstring)
    - juli 21: add phred as list into seq_list in plot_base_grid
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec



def plot_base_grid(seq_list: list,
                   x_labels: list = None,
                   y_labels: list = None,
                   height: float = 4,
                   width: float = 8,
                   colors: dict = None,
                   fig_dir: str = './figure/',
                   fig_name: str = 'viz_bases',
                   plot_title: str = "DNA sequence",
                   x_lab: str = "Position",
                   y_lab: str = "Sequence",
                   fsize: float = 12,
                   artificial=False
                   ):
    """
    plot the headmap-style of list of DNA sequences
    
    each sequence is displayed as a row in a grid, with each nucleotide represented
    by a colored cell labeled with the base character. 
    
    Parameters
    ----------
    seq_list : 
        A list of DNA sequences to plot. Each sequence should be the same length.
        
    x_labels : 
        Labels for the x-axis (positions). Defaults to numeric indices.
        
    y_labels : 
        Labels for the y-axis (sequence IDs). Defaults to numeric indices.
        
    height : height of figure
        
    width : width of figure
        
    colors : 
        dictionary mapping base characters (A, T, C, G, N, -) to color hex codes. 
        
    fig_dir : 
        directory where the figure will be saved. 
        
    fig_name :
        name of the saved figure file

    plot_title :
        title for the plot
        
    x_lab :
        label for the x-axis
        
    y_lab : str, optional
        label for the y-axis
    
    fsize : size of text

    Returns
    -------
    None
        The function saves the plot as an image file and does not return any value.

    Raises
    ------
    IndexError
        If the sequences in `seq_list` are of different lengths.

    Examples
    --------
    >>> seq_list = [
            'GACAAACCAGAAGCCGCTCC',
            'GTTCACACCCATGACGAACA',
            'GAACACAAAGCATAGACTGC',
            'GAGTCCGAGCAGAAGAAGAA',
            'GAAGACCAAGGATAGACTGC',
            'GAA-ACCAAGGATAGACTNC'
            ]
    >>> plot_base_grid(seq_list)
    """
    
    # Validate equal length of all sequences
    seq_lengths = [len(seq) for seq in seq_list]
    if len(set(seq_lengths)) > 1:
        if artificial:
            max_len = max(seq_lengths)
            seq_list[:] = [
                seq + 'X' * (max_len - len(seq)) if isinstance(seq, str) else 
                seq + ['X'] * (max_len - len(seq)) if isinstance(seq, list) else 
                seq
                for seq in seq_list]
        else:
            raise ValueError("All sequences in `seq_list` must be of the same length.")
    # convert sequences to uppercase character lists
    #seq_list_bases = [list(seq.upper()) for seq in seq_list]
    seq_list_bases = [
        list(seq.upper()) if isinstance(seq, str) else 
        [str(base).upper() for base in seq] if isinstance(seq, list) else 
        [str(base) for base in seq]
        for seq in seq_list
    ]
    
    if colors is None:
        colors = {'A': '#FF9999', 'a': '#FF6666',
                  'T': '#99FF99', 't': '#66FF66',
                  'C': '#9999FF', 'c': '#6666FF',
                  'G': '#FFFF99', 'g': '#FFFF66',
                  '-': '#FFFFFF',
                  'N': '#A0A0A0', 'n': '#808080',
                  '×': '#1C2526', 'X': '#0E1213'
                  }
        for i in range(1, 41):
            t = (i - 1) / 39.0
            # White (#FFFFFF, RGB: 255,255,255) to orange (#FF4500, RGB: 255,69,0)
            r = 255  # Red stays constant
            g = int(255 - (255 - 69) * t)  # 255 to 69
            b = int(255 - (255 - 0) * t)   # 255 to 0
            colors[str(i)] = f'#{r:02X}{g:02X}{b:02X}'
    
    # Create figure and axis
    fig = plt.figure(figsize=(width, height))
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0, 0])
     
    # Plot bases with colored rectangles
    for i, row in enumerate(seq_list_bases):
        for j, base in enumerate(row):
            rect = plt.Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor=colors.get(base, '#CCCCCC'))
            ax.add_patch(rect)
            # text color: white for PHRED scores (strings of integers 1-40), black otherwise
            text_color = 'white' if base.isdigit() and 1 <= int(base) <= 40 else 'black'
            #ax.text(j, i, base, ha='center', va='center', color='black', fontsize=fsize)
            ax.text(j, i, base, ha='center', va='center', color=text_color, fontsize=fsize)
    # #
    # num_base = len(seq_list_bases[0])
    # num_seq = len(seq_list_bases)
    # # Axis labels and ticks
    # if x_labels is None:
    #     x_labels = list(range(num_base))
    
    # if y_labels is None:
    #     y_labels = list(range(num_seq))
    
    # x_axis_labels = [label if label in x_labels else '' for label in x_labels]
    
    # # Set limits and grid
    # ax.set_aspect('equal', adjustable='box')
    # ax.set_xlim(-0.5, num_base - 0.5)
    # ax.set_ylim(-0.5, num_seq - 0.5)
    # ax.set_xticks(range(len(x_labels)))
    # ax.set_xticklabels(x_axis_labels)
    # #ax.set_yticks(range(len(y_labels)))
    # ax.set_yticklabels(y_labels)
    # ax.set_title(plot_title)
    # ax.set_xlabel(x_lab)
    # ax.set_ylabel(y_lab)
    # ax.invert_yaxis()
    # plt.tight_layout()
    # Get dimensions
    num_base = len(seq_list_bases[0])
    num_seq = len(seq_list_bases)
    
    # set x-axis labels
    if x_labels is None:
        x_labels = list(range(num_base))    # default to 0-based indices
    ax.set_xticks(range(num_base))          # one tick per base position
    x_axis_labels = [pos if pos in x_labels else '' for pos in range(num_base)]
    ax.set_xticklabels(x_axis_labels, rotation=45, ha='right')       # assign labels directly
    
    # set y-axis labels
    if y_labels is None:
        y_labels = list(range(num_seq))
    ax.set_yticks(range(num_seq))
    ax.set_yticklabels(y_labels)
    
    # set plot properties
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(-0.5, num_base - 0.5)
    ax.set_ylim(-0.5, num_seq - 0.5)
    ax.set_title(plot_title)
    ax.set_xlabel(x_lab)
    ax.set_ylabel(y_lab)
    ax.invert_yaxis()
    plt.tight_layout()
    
    # save figure
    os.makedirs(fig_dir, exist_ok=True)
    plt.savefig(os.path.join(fig_dir, fig_name))
    plt.close()


def cigar_recover(sequence, cigar):
    """
    Recover the original sequence by interpreting the CIGAR string.
    
    CIGAR operations are handled as follows:
    - M: Keep the sequence as is.
    - D: Add '-' for deleted bases.
    - I: Convert bases to lowercase.
    - S: Keep the sequence as is (soft-clipped bases are retained).
    - H: Add '×' for hard-clipped bases.
    
    Args:
        sequence (str): The input sequence (e.g., 'CCTAGTCCAAACTGGATCTCTGCTGTCCCTG').
        cigar (str): The CIGAR string (e.g., '224H31M').
    
    Returns:
        str: The recovered sequence with modifications based on CIGAR operations.
    
    Examples:
        >>> sequence = 'CCTAGTCCAAACTGG'
        >>> cigar = '5M'  # Keep 5 bases as is
        >>> cigar_recover(sequence, cigar)
        'CCTAG'
        
        >>> cigar = '3H5M'  # Add 3 '×' then keep 5 bases
        >>> cigar_recover(sequence, cigar)
        '×××CCTAG'
        
        >>> cigar = '5M3D'  # Keep 5 bases, add 3 '-'
        >>> cigar_recover(sequence, cigar)
        'CCTAG---'
        
        >>> cigar = '5M3I'  # Keep 5 bases, lowercase 3 bases
        >>> cigar_recover(sequence, cigar)
        'CCTAGtcc'
        
        >>> cigar = '5M3S'  # Keep 5 bases, keep 3 bases (soft-clipped)
        >>> cigar_recover(sequence, cigar)
        'CCTAGTCC'
        
        >>> cigar = '2H3M2D2I3S'  # Combined: 2 '×', 3 bases, 2 '-', 2 lowercase, 3 bases
        >>> cigar_recover(sequence, cigar)
        '××CCT--aaACT'
        
        >>> sequence = 'CCTAGTCCAAACTGGATCTCTGCTGTCCCTG'
        >>> cigar = '224H31M'
        >>> cigar_recover(sequence, cigar)
        '××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××CCTAGTCCAAACTGGATCTCTGCTGTCCCTG'
    """
    import re
    
    # Parse CIGAR string into operations (e.g., [('224', 'H'), ('31', 'M')])
    cigar_ops = re.findall(r'(\d+)([HMSID])', cigar)
    
    result = []
    seq_pos = 0  # Position in the input sequence
    
    for count, op in cigar_ops:
        count = int(count)
        if op == 'H':
            # Add 'count' number of '×' for hard-clipped bases
            result.append('×' * count)
        elif op == 'D':
            # Add 'count' number of '-' for deleted bases
            result.append('-' * count)
        elif op == 'M' or op == 'S':
            # Keep 'count' bases from the sequence as is
            result.append(sequence[seq_pos:seq_pos + count])
            seq_pos += count
        elif op == 'I':
            # Convert 'count' bases to lowercase
            result.append(sequence[seq_pos:seq_pos + count].lower())
            seq_pos += count
    
    return ''.join(result)
