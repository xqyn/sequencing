'''
project: sequencing
april 29 2025 - XQ - Leiden UMC
making barplot
''' 
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import LogLocator
import numpy as np
import pandas as pd


# --- plot_bar
def plot_bar(df,
             title='add_title',             
             color_bar=None,
             height=7,
             width=12,
             bar_width=0.25,
             xlabel='samples',
             ylabel='log_count',
             ysticks=None,
             gap=None,
             fig_dir=None
             ):
    """
    Create a bar plot for unique barcode pairing spacers from a DataFrame.
    
    Parameters
    ----------
    df : pandas.DataFrame
        dataFrame with days as index and series names as columns.
    title : str
        title of the plot (default: 'add_title').
    color_bar : list
        List of colors for bars (default: None, uses tableau colors).
    height : float
        height of the figure in inches (default: 7).
    width : float
        width of the figure in inches (default: 12).
    bar_width : float
        width of the bars in the plot (default: 0.25).
    xlabel : str
        label for the x-axis (default: 'samples').
    ylabel : str
        label for the y-axis (default: 'log_count').
    ysticks : list
        custom y-axis tick labels (default: None).
    gap : float
        gap between bars and text value (default: None).
    fig_dir : str
        directory path to save the figure (default: None).
    
    Returns
    -------
    matplotlib.figure.Figure
        The generated figure object for further customization if needed.
    
    Example
    -------
    >>> bc_sp_pair = pd.DataFrame({
    ...     'ABE': [21431, 12758, 17968, 20056, 11116, 15565],
    ...     'CBE': [21074, 27101, 21821, 31957, 16373, 23149]
    ... }, index=['d0', 'd1', 'd3', 'd5', 'd10', 'd10u'])
    >>> ystick = [2.5e3, 5e3, 7.5e3, 1e4, 2.5e4, 5e4]
    >>> fig_dir = './figures/'
    >>> fig = plot_bar(bc_sp_pair, title='barcode_spacer_pair', 
    ...                color_bar=['#f7b5a8', '#f06449'], 
    ...                fig_dir=fig_dir, height=8, width=12, 
    ...                ysticks=ystick)
    """
    
    # input validation:
    if not isinstance(df, pd.DataFrame) or df.empty:
        raise ValueError('input "df" must be non-empty panda DataFrame')
    
    # default colors if none provided
    color_bar = plt.get_cmap('tab10').colors if color_bar is None else color_bar
    
    # close any existing plots
    plt.close()
    
    # create base figure and axs
    fig, ax = plt.subplots(figsize=(width, height))
    
    # bar position
    bar_pos = np.arange(len(df.index))
    n_series = len(df.columns)
    #bar_width = 0.8 / n_series  # Adjust width based on number of series
    
    # # gap for text above the bar plot
    # gap = df.mean().mean() * 0.05 if gap is None else gap
    # plot bar for each series
    for idx, series in enumerate(df.columns):
        values = df[series].values
        offset = (idx - n_series / 2 + 0.5) * bar_width
        positions = bar_pos + offset
        
        # plot bars
        ax.bar(
            positions, 
            values, 
            width=bar_width,
            label=series, 
            color=color_bar[idx % len(color_bar)])
        
        # gap for text above the bar plot, specific to each column
        gap_c = df[series].mean() * 0.05 if gap is None else gap
        # Annotate bars with count values
        for pos, value in zip(positions, values):
            ax.text(pos, 
                    value + gap_c, 
                    str(value),
                    ha='center', 
                    va='bottom', 
                    rotation=90, 
                    fontsize=10)
        
    # labels and title
    ax.set_title(title)
    ax.legend(loc='upper right', bbox_to_anchor=(1.1, 1))
    
    # x
    ax.set_xlabel(xlabel)
    ax.set_xticks(bar_pos)
    ax.set_xticklabels(df.index, rotation=45)
    
    # y
    ax.set_ylabel(ylabel)
    ax.set_yscale("log")
    # dynamic y-ticks if not provided
    if ysticks is None:
        max_value = df.max().max()
        min_value = df[df > 0].min().min() if (df > 0).any().any() else 1
        log_min = int(np.floor(np.log10(min_value)))
        log_max = int(np.ceil(np.log10(max_value)))
        
        # Clean powers of 10
        ysticks = np.logspace(log_min - 0.5, log_max, log_max - log_min + 2)
    
    ax.set_yticks(ysticks)
    ax.set_yticklabels([f"${y/10**int(np.log10(y)):.1f} \\times 10^{{{int(np.log10(y))}}}$" if y != 10**int(np.log10(y)) else f"$10^{{{int(np.log10(y))}}}$" for y in ysticks])
    #ax.yaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,), numticks=10))
    
    # grid
    ax.grid(True, linestyle='--', alpha=0.5) 
    # Enable minor ticks
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs='auto', numticks=100))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())  # hide minor labels
    ax.set_axisbelow(True)
    # Gridlines
    ax.grid(True, which='major', linestyle='--', alpha=0.75)
    ax.grid(True, which='minor', linestyle=':', alpha=0.5)
    
    # remove top and right spines (from original function)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # save figure
    if fig_dir is not None:
        plt.savefig(f"{fig_dir}{title}.png", bbox_inches='tight');
    plt.close()  # Close the figure
    # return figure for further customization if needed
    return fig