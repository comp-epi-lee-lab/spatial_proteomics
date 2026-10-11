from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from spatial_proteomics.cell_proportions import _proportion_by_folder_name

import logging
logger = logging.getLogger(__name__)


def _pct_dists(len_sizes):
    """
    Calculate radial positions for pie plot percentage labels.

    Parameters
    ----------
    len_sizes : int
        Number of pie plot slices. Intended for values greater than one.

    Returns
    -------
    list of float
        Radial scaling factors used to separate percentage labels.
    """
    splitting_number = (1.05-(1/len_sizes))/(len_sizes-1)
    distances_list = np.round(np.arange(1/len_sizes, 1.05+.01,splitting_number),2).tolist()
    just_because = len(distances_list)
    return [distances_list[just_because-1-idx] for idx in range(just_because)]

def _move_pct(t, d):
    """
    Move a Matplotlib pie plot percentage label radially.

    Parameters
    ----------
    t : matplotlib.text.Text
        Percentage label created by ``Axes.pie``.
    d : float
        Radial scaling factor.

    Returns
    -------
    None
        The text object is modified in place.
    """
    xi, yi = t.get_position()
    ri = np.sqrt(xi**2 + yi**2)
    phi = np.arctan2(yi, xi)
    x = d * ri * np.cos(phi)
    y = d * ri * np.sin(phi)
    t.set_position((x, y))

def _plot_pie_plots(sum_dict,custom_colors,title_name,unique_values,ct_name,save_namefile):
    """
    Save a pie plot showing cell population composition for one sample.

    Parameters
    ----------
    sum_dict : dict
        Cell counts keyed by phenotype; must include ``"Total cells"``.
    custom_colors : dict of str to str
        Phenotype-to-color mapping.
    title_name : str
        Figure title.
    unique_values : list of str
        Populations represented in the sample.
    ct_name : str
        Focal population placed first in the pie plot.
    save_namefile : pathlib.Path
        Output image path.

    Returns
    -------
    None
        The figure is written below the ``Statistics`` output directory.
    """
    fig, ax = plt.subplots(figsize=(8, 8), constrained_layout=True)
    
    k1=ct_name
    if sum_dict['Total cells'] > 0 and sum_dict[k1] > 0:                # first cell type will be the "mouth" of this 'pac-man' plot if it is 'open'
        blue_frac = sum_dict[k1] / sum_dict['Total cells']
        half_blue_angle_deg = (blue_frac * 360) / 2.0
        startangle = half_blue_angle_deg
    else:
        startangle = 90
        
    colors = [custom_colors[k1]]+[color for k,color in custom_colors.items() if k in unique_values and k != k1]
    sizes  = [sum_dict[k1]]+[sum_dict[cell_type] for cell_type,_ in custom_colors.items() if cell_type in unique_values and cell_type != k1]
    pctdists1 = _pct_dists(len(sizes)) if len(sizes) > 1 else [1]                                   # here, we are equally distributing labels across slices of pie plot

    wedges, texts, autotexts = ax.pie(
        sizes,
        colors=colors,
        startangle=startangle,
        counterclock=False,
        autopct=(lambda p: f'{p:.2f}%') if sum_dict['Total cells'] else None,
        pctdistance=0.75,
        wedgeprops=dict(linewidth=1, edgecolor='white')
    )

    ax.set_title(title_name)
    for t in texts: t.set_visible(False)
    for t, d in zip(autotexts, pctdists1): _move_pct(t, d)              # here, we are putting labels using pct dists
    ax.axis('equal')

    ax.legend(
        handles=wedges,
        labels = [f"{k1} ({sum_dict[k1]})"]+[f"{cell_type} ({sum_dict[cell_type]})" for cell_type,_ in custom_colors.items() if cell_type in unique_values and cell_type != k1],
        loc='lower center', bbox_to_anchor=(0.5, -0.15), ncol=1, frameon=False
    )
    
    plt.savefig(save_namefile)
    plt.close()


def calculate_cell_proportion_matrix(adata_dicts, output_dir, custom_colors, focus_ct_name_list, pie_plot=False,overwrite_existing_files=False):
    """
    Build sample-by-phenotype proportion matrices and optional pie plots.

    Parameters
    ----------
    adata_dicts : dict of str to anndata.AnnData
        Per-sample annotated data.
    output_dir : str or pathlib.Path
        SpPrAn output directory.
    custom_colors : dict of str to str
        Phenotype-to-color mapping.
    focus_ct_name_list : list of str
        Major cell types prioritized as focal populations in pie plots.
    pie_plot : bool, default=False
        Whether sample-level composition pie plots should be generated.
    overwrite_existing_files : bool, default=False
        Whether existing outputs may be replaced.

    Returns
    -------
    None
        Proportion matrices and optional pie plots are written below the
        ``Statistics`` output directory.
    """
    output_dir = Path(output_dir)
    results_dir = output_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    
    prop_matrix_list = {}
    
    for k, adata in adata_dicts.items():
        if adata.n_obs == 0:
            logger.warning(f"No cells were found in Sample {k}. Skipping sample...")
            continue
        if "color_keys" in adata.uns and "color_vals" in adata.uns:
            for key,val in zip(adata.uns['color_keys'], adata.uns['color_vals']):
                if key not in custom_colors.keys(): custom_colors[key] = val

        for folder_name in adata.uns['folder_names']:
            if folder_name not in prop_matrix_list.keys(): prop_matrix_list[folder_name]=[]
            _, df, title_name = _proportion_by_folder_name(custom_colors, overwrite_existing_files, results_dir, k, adata, folder_name)
            if df is None: continue
            df.index = ["Total cells in cell type", "Percentage"]
            prop_matrix_list[folder_name].append([k,df])
            if pie_plot:
                statistics_plots_dir = output_dir / "Statistics" / "plots" / folder_name / "pie_plots"
                statistics_plots_dir.mkdir(parents=True, exist_ok=True)
                pie_title_name = f"{title_name} - {folder_name} - Pie plot.png"
                save_namefile = statistics_plots_dir / pie_title_name
                if save_namefile.exists() and not overwrite_existing_files: 
                    logger.warning(f"File '{pie_title_name}' already exists...")
                    continue
                sum_dict = {col_name:int(df[col_name]["Total cells in cell type"]) for col_name in df.columns}
                unique_values = df.columns.tolist()
                ct_name = unique_values[0]
                for ctname in focus_ct_name_list:
                    if ctname in df.columns.tolist():
                        ct_name = ctname
                        break
                sum_dict['Total cells'] = int(df.loc["Total cells in cell type"].sum())
                _plot_pie_plots(sum_dict, custom_colors, title_name = pie_title_name, unique_values = unique_values, ct_name = ct_name, save_namefile=save_namefile)
    for key, val_list in prop_matrix_list.items():
        statistics_results_dir = output_dir / "Statistics" / "results" /  key 
        statistics_results_dir.mkdir(parents=True, exist_ok=True)
        csv_filename = statistics_results_dir / "Cell type proportions matrix.csv"
        if csv_filename.exists() and not overwrite_existing_files:
            logger.warning(f"File 'Cell type proportions matrix.csv' for {key} already exists...")
            continue
        for i, df_to_concat in enumerate(val_list):
            s = df_to_concat[1].loc["Percentage"].copy()
            s.name = df_to_concat[0]   # sample id
            prop_matrix = s if i==0 else pd.concat([prop_matrix,s],axis=1)
        prop_matrix = prop_matrix.T
        prop_matrix.fillna("0%", inplace=True)
        prop_matrix.to_csv(csv_filename)
    if pie_plot: logger.warning("Pie plots and cell type proportion matrix generated!")
    else: logger.warning("Cell type proportion matrix generated!")