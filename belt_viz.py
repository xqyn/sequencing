"""
project: belt
XQ - Leiden UMC
visualization of squencing alignments
update: 2025
    - juni 16: plot_base_grid
    - juni 26: fix x_label for position base + fsize adds
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
                   fsize: float = 12
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
        raise ValueError("All sequences in `seq_list` must be of the same length.")
    
    # convert sequences to uppercase character lists
    seq_list_bases = [list(seq.upper()) for seq in seq_list]
    
    
    if colors is None:
        colors = {'A': '#FF9999', 
                  'T': '#99FF99', 
                  'C': '#9999FF', 
                  'G': '#FFFF99',
                  '-': '#FFFFFF',
                  'N': '#808080',
                  }
    
    # Create figure and axis
    fig = plt.figure(figsize=(width, height))
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0, 0])
    
    # Plot bases with colored rectangles
    for i, row in enumerate(seq_list_bases):
        for j, base in enumerate(row):
            rect = plt.Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor=colors.get(base, '#CCCCCC'))
            ax.add_patch(rect)
            ax.text(j, i, base, ha='center', va='center', color='black', fontsize=fsize)
     
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
    
    # Set x-axis labels
    if x_labels is None:
        x_labels = list(range(num_base))  # Default to 0-based indices
    
    # # Ensure x_labels length matches num_base
    # if len(x_labels) != num_base:
    #     x_labels = list(range(num_base))  # Fallback to default if mismatch
    
    # Set ticks and labels
    ax.set_xticks(range(num_base))  # One tick per base position
    x_axis_labels = [pos if pos in x_labels else '' for pos in range(num_base)]
    ax.set_xticklabels(x_axis_labels)    # Assign labels directly
    
    # Set y-axis labels
    if y_labels is None:
        y_labels = list(range(num_seq))  # Default to 0-based indices
    
    ax.set_yticks(range(num_seq))    # One tick per sequence
    ax.set_yticklabels(y_labels)     # Assign labels directly
    
    # Set plot properties
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(-0.5, num_base - 0.5)
    ax.set_ylim(-0.5, num_seq - 0.5)
    ax.set_title(plot_title)
    ax.set_xlabel(x_lab)
    ax.set_ylabel(y_lab)
    ax.invert_yaxis()
    plt.tight_layout()
    
    # Save figure
    os.makedirs(fig_dir, exist_ok=True)
    plt.savefig(os.path.join(fig_dir, fig_name))
    plt.close()




