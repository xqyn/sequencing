'''
project: sequencing
april 29 2025 - XQ - Leiden UMC
making barplot
''' 
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


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
    
    Parameters:
    df (pandas.DataFrame): DataFrame with days as index and series names as columns
    fig_dir (str): Directory path to save the figure
    colour_bar (list, optional): List of colors for bars. Defaults to tableau colors
    """
    # input validation:
    if not isinstance(df, pd.DataFrame) or df.empty:
        raise ValueError('input "df" must be non-empty panda DataFrame')
    
    # default colors if none provided
    color_bar = plt.get_cmap('tab10').colors if color_bar is None else color_bar
    
    # gap for text above the bar plot
    gap = df.mean().mean() * 0.05 if gap is None else gap
    
    # close any existing plots
    plt.close()
    
    # create base figure and axs
    fig, ax = plt.subplots(figsize=(width, height))
    
    # bar position
    bar_pos = np.arange(len(df.index))
    n_series = len(df.columns)
    #bar_width = 0.8 / n_series  # Adjust width based on number of series
    
    # plot bar for each series
    for idx, series in enumerate(df.columns):
        values = df[series].values
        offset = (idx - n_series / 2 + 0.5) * bar_width
        positions = bar_pos + offset
        
        # Ppot bars
        ax.bar(
            positions, 
            values, 
            width=bar_width,
            label=series, 
            color=color_bar[idx % len(color_bar)])
        
        # Annotate bars with count values
        for pos, value in zip(positions, values):
            ax.text(pos, 
                    value + gap, 
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
    ax.set_xticklabels(df.index)
    
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
    ax.set_yticklabels([f"$10^{int(np.log10(y))}$" for y in ysticks])
    
    # grid
    ax.grid(True, linestyle='--', alpha=0.5) 
    # Enable minor ticks
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs='auto', numticks=100))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())  # hide minor labels
    ax.set_axisbelow(True)
    # Gridlines
    ax.grid(True, which='major', linestyle='--', alpha=0.6)
    ax.grid(True, which='minor', linestyle=':', alpha=0.3)
    
    # remove top and right spines (from original function)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # save figure
    if fig_dir is not None:
        plt.savefig(f"{fig_dir}{title}.png",
                    bbox_inches='tight')
    
    # return figure for further customization if needed
    return fig