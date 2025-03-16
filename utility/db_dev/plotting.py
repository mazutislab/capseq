from anndata import AnnData
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def condition_summary(data, groupby=[], output_type="absolute"):
    """
    Summarizes cell counts or percentages based on specified groupings.
    
    Parameters:
        data (AnnData or pd.DataFrame): Input data, either an AnnData object or a DataFrame.
        groupby (list): List of column names to group by.
        output_type (str): "absolute" for absolute counts, "percentage" for row-wise percentages.
        
    Returns:
        pd.DataFrame: DataFrame containing either absolute counts or row-wise percentages.
    """

    if isinstance(data, AnnData):
        obs_data = data.obs
    elif isinstance(data, pd.DataFrame):
        obs_data = data
    else:
        raise TypeError("Input data must be an AnnData object or a pandas DataFrame.")
    

    if not all(col in obs_data.columns for col in groupby):
        raise ValueError("All columns in groupby must be present in the data.")
    
    if len(groupby) == 1:
        cell_counts = pd.DataFrame(obs_data.groupby(groupby).size(), columns=["cell_number"])
    else:
        cell_counts = obs_data.groupby(groupby).size().unstack()
    
    if output_type == "percentage":
        if len(groupby) == 1:
            raise ValueError("Row-wise percentages require at least two grouping columns.")
        else:
            # Normalize each row to get row-wise percentages
            cell_counts = cell_counts.div(cell_counts.sum(axis=1), axis=0) * 100

    return cell_counts


def hbars(results, category_names, category_colors, figsize=(5, 3), plot_legend = True):
    """
    Parameters
    ----------
    results : dict
        A mapping from question labels to a list of answers per category.
        It is assumed all lists contain the same number of entries and that
        it matches the length of *category_names*.
    category_names : list of str
        The category labels.
    """
    labels = list(results.keys())
    data = np.array(list(results.values()))
    data_cum = data.cumsum(axis=1)
    
    fig, ax = plt.subplots(figsize=figsize, dpi=300)
    
    ax.invert_yaxis()
    ax.xaxis.set_visible(False)
    ax.set_xlim(0, np.sum(data, axis=1).max())

    # Removing the frame of the figure
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)

    ax.grid(False)
    
    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        widths = data[:, i]
        starts = data_cum[:, i] - widths
        ax.barh(labels, widths, left=starts, height=0.9, label=colname, color=color)
        xcenters = starts + widths / 2
    
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(6)

    if plot_legend == True:
        ax.legend(
            loc='center left',
            bbox_to_anchor=(1, 0.5),
            fontsize = 6,
            frameon = False
        )

    plt.tight_layout()
    return fig, ax