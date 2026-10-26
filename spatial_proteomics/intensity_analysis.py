import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging

logger = logging.getLogger(__name__)

def _read_intensities(adata_dicts, which_pm = 'all'):
    """
    Extract protein marker intensity matrices from AnnData objects.

    Parameters
    ----------
    adata_dicts : dict of str to anndata.AnnData
        Per-sample AnnData objects.
    which_pm : "all" or list of str, default="all"
        AnnData variable names to retain.

    Returns
    -------
    dict of str to pandas.DataFrame
        Per-sample intensity matrices indexed by cell identifier.
    """
    intensities_dict = {}
    for k, adata in adata_dicts.items():
        if which_pm == 'all': intensities_dict[k] = pd.DataFrame(adata.X,columns=adata.var_names,index=adata.obs_names)
        else: intensities_dict[k] = pd.DataFrame(adata.X,columns=adata.var_names,index=adata.obs_names)[which_pm]
    return intensities_dict

def _boxplot_calculation(df,pm_name):
    """
    Calculate Tukey box plot summary values for one protein marker.

    Returns
    -------
    list of float
        Lower whisker, Q1, median, Q3, and upper whisker using the 1.5-IQR rule.
    """
    q1 = df[pm_name].quantile(0.25)
    median = df[pm_name].quantile(0.50)
    q3 = df[pm_name].quantile(0.75)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    inner_data = df[(df[pm_name] >= lower_bound) & (df[pm_name] <= upper_bound)][pm_name]
    lower_whisker = inner_data.min()
    upper_whisker = inner_data.max()
    return [lower_whisker, q1, median, q3, upper_whisker]

def _boxplot_outliers(df,pm_name):
    """
    Return intensity observations outside the 1.5-IQR whisker range.

    Returns
    -------
    pandas.Series
        Sorted outlier intensity values.

    Notes
    -----
    It is currently not used by ``_intensities_data``.
    """
    q1 = df[pm_name].quantile(0.25)
    q3 = df[pm_name].quantile(0.75)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    outliers = df[(df[pm_name] < lower_bound) | (df[pm_name] > upper_bound)][pm_name].sort_values()
    return pd.Series(outliers)

def _intensities_data(intensities_dict,adata_dicts,list_ct_short_name,output_dir):
    """
    Export box plot summary data by discovered phenotype.

    Parameters
    ----------
    intensities_dict : dict of str to pandas.DataFrame
        Per-sample protein marker intensity matrices.
    adata_dicts : dict of str to anndata.AnnData
        AnnData objects containing discovered phenotype labels.
    list_ct_short_name : list of str
        Short names of explored parent populations.
    output_dir : pathlib.Path
        Base directory for intensity analysis results.

    Returns
    -------
    None
        Per-sample/per-phenotype Tukey box-plot summaries and complete intensity-data 
        tables are written as CSV files.

    Notes
    -----
    Long marker-defined phenotype names are mapped to compact IDs such as
    ``subTumor_P001`` for output filenames. A per-sample phenotype-mapping CSV
    preserves the correspondence between compact IDs and complete phenotype names.
    """

    for k, adata in adata_dicts.items():
        intensities_csv_dir_phen = output_dir / "results" / "Intensity analysis by phenotype" / f"{k}"
        intensities_csv_dir_phen.mkdir(parents=True, exist_ok=True)
        intensity_df = intensities_dict[k]
        pm_list = intensity_df.columns.tolist()

        for short_ct_name in list_ct_short_name:
            host_column = f"sub{short_ct_name} phenotype"
            
            list_host_names = adata.obs[host_column].unique().tolist()
            n_digits = max(3, len(str(len(list_host_names))))
            phenotype_ids = {phenotype: f"sub{short_ct_name}_P{i:0{n_digits}d}" for i, phenotype in enumerate(list_host_names)}
            mapping = pd.DataFrame({"Phenotype ID": [pvals for pvals in phenotype_ids.values()], "Phenotype": [pkeys for pkeys in phenotype_ids.keys()]})
            intensities_csv_dir_phen_ext = intensities_csv_dir_phen / short_ct_name
            intensities_csv_dir_phen_ext.mkdir(parents=True, exist_ok=True)
            mapping.to_csv(intensities_csv_dir_phen / f"Phenotype mapping - Sample {k} - sub{short_ct_name}.csv", index=False)
            mapping = mapping.set_index('Phenotype')
            for phenotype in list_host_names:
                if pd.isna(phenotype):
                    continue

                cell_ids = adata.obs_names[adata.obs[host_column] == phenotype].tolist()
                filtered_df = intensity_df[intensity_df.index.isin(cell_ids)]
                text_total_cells = "(1 cell)" if len(filtered_df)==1 else f"({len(filtered_df)} cells)"

                title_name_1 = f"Box plot data - Sample {k} - {host_column} - {mapping.loc[phenotype,"Phenotype ID"]} {text_total_cells}"
                save_namefile_1 = intensities_csv_dir_phen_ext / f"{title_name_1}.csv"
                # title_name_2 = f"Box plot data - Sample {k} - {host_column} - {phenotype} - outliers"
                # save_namefile_2 = intensities_csv_dir_phen_ext / f"{title_name_2}.csv"

                title_name_3 = f"Intensity data - Sample {k} - {host_column} - {mapping.loc[phenotype,"Phenotype ID"]} {text_total_cells}"
                save_namefile_3 = intensities_csv_dir_phen_ext / f"{title_name_3}.csv"

                pd.DataFrame({'Metric': ['Lower_Whisker', 'Q1', 'Median (Q2)', 'Q3', 'Upper_Whisker']}|{ pm_name: _boxplot_calculation(filtered_df,pm_name) for pm_name in pm_list }).to_csv(save_namefile_1, index=False)
                # pd.DataFrame({ pm_name: _boxplot_outliers(filtered_df,pm_name) for pm_name in pm_list }).to_csv(save_namefile_2, index=False)
                pd.DataFrame({ pm_name: filtered_df[pm_name] for pm_name in pm_list }).to_csv(save_namefile_3, index=False)

def calculate_pm_intensity_analysis(adata_dicts, list_ct_short_name, output_dir, which_pm):
    """
    Export protein marker intensity summaries for discovered phenotypes.

    Parameters
    ----------
    adata_dicts : dict of str to anndata.AnnData
        AnnData objects containing discovered phenotype annotations.
    list_ct_short_name : list of str
        Short names of explored parent populations.
    output_dir : pathlib.Path
        SpPrAn output directory.
    which_pm : list of str
        Positivity column names whose corresponding intensity measurements
        should be summarized. An empty list disables the analysis.

    Returns
    -------
    None
        Summary tables are written under the statistics/intensity output tree.

    Notes
    -----
    This is descriptive output for box plot visualization; the function itself
    does not perform inferential tests on marker intensity.
    """
    if len(which_pm)!=0:
        positivity_column_names = pd.read_csv(output_dir / "results" / "Positivity column names.csv")
        epm_list = positivity_column_names[positivity_column_names["Positivity column names"].isin(which_pm)]["Name appearing in anndata file"].tolist()
        intensities_dict = _read_intensities(adata_dicts, epm_list)
        _intensities_data(intensities_dict, adata_dicts, list_ct_short_name, output_dir / "Statistics")

def _getting_vmin_n_vmax (adata_dicts):
    """
    Calculate cross-sample protein-marker-intensity limits from stored thresholds.

    Parameters
    ----------
    adata_dicts : dict of str to anndata.AnnData
        Samples containing protein marker thresholds in ``adata.uns["thresholds"]``.

    Returns
    -------
    dict of str to tuple of float
        Global minimum/maximum threshold limits for each protein marker.

    Notes
    -----
    This helper is intended to support consistent color scales across samples.
    It is currently not used by ``plot_spatial_protein_marker_intensities``.
    """
    marker_limits = {}
    positivity_columns = []
    for k, adata in adata_dicts.items():
        # adata.uns["thresholds"] ={ "pm_intensities" : { "low_treshold": "vaue", high_treshold" : "value"} }
        if isinstance(adata.uns["thresholds"], str): continue
        if len(positivity_columns) == 0:
            positivity_columns = list(adata.var_names)
            for pm in positivity_columns:
                marker_limits[pm] = (adata.uns["thresholds"][pm]["low_treshold"], adata.uns["thresholds"][pm]["high_treshold"])
        else:
            for pm in positivity_columns:
                vmin_temp, vmax_temp = marker_limits[pm]
                marker_limits[pm] = (
                    min(vmin_temp, adata.uns["thresholds"][pm]["low_treshold"]),
                    max(vmax_temp, adata.uns["thresholds"][pm]["high_treshold"])
                )
    return marker_limits

def plot_spatial_protein_marker_intensities(adata_dicts, output_dir, cmap="viridis", point_size=1, dpi=150):
    """
    Plot spatial protein-marker intensity for every sample and protein marker.

    Parameters
    ----------
    adata_dicts : dict of str to anndata.AnnData
        Samples containing spatial coordinates in ``obsm["spatial"]`` and
        protein-marker intensities in ``X``.
    output_dir : pathlib.Path
        SpPrAn output directory.
    cmap : str, default="viridis"
        Matplotlib colormap used for intensity values.
    point_size : float, default=1
        Scatter-point size for individual cells.
    dpi : int or float, default=150
        Output image resolution.

    Returns
    -------
    None
        One PNG per sample and protein marker is written under
        ``plots/Spatial protein markers intensity plots``.

    Notes
    -----
    Each plot colors individual segmented cells by continuous protein-marker
    intensity, analogous to displaying gene expression on a single-cell
    embedding but retaining the cells' actual tissue coordinates.

    The current implementation scales each sample/protein marker independently using
    its observed minimum and maximum intensity. Consequently, colors should not
    be compared quantitatively between different sample panels unless a common
    scale is introduced.
    """
    logger.warning(f"Generating protein markers' intensities spatial plots...")
    intensity_dir = output_dir / "plots" / "Spatial protein markers intensity plots"
    intensity_dir.mkdir(parents=True, exist_ok=True)
    # marker_limits = _getting_vmin_n_vmax (adata_dicts)
    marker_limits = {}
    for k, adata in adata_dicts.items():
        intensity_dir_2 = intensity_dir / f"{k}"
        intensity_dir_2.mkdir(parents=True, exist_ok=True)
        coords = adata.obsm["spatial"]
        matrix = adata.X
        for idx, marker in enumerate(adata.var_names):
            intensities = np.asarray(matrix[:, idx]).ravel()
            if len(marker_limits)!=0 and marker in marker_limits: vmin, vmax = marker_limits[marker]
            else:
                vmin = np.nanmin(intensities)
                vmax = np.nanmax(intensities)
            fig, ax = plt.subplots(figsize=(6, 6))
            scatter = ax.scatter(coords[:, 0], coords[:, 1], c=intensities, cmap=cmap, s=point_size, vmin=vmin, vmax=vmax, edgecolors="none", linewidths=0, rasterized=True)
            ax.set_aspect("equal")
            ax.invert_yaxis()
            ax.axis("off")
            ax.set_title(marker)
            fig.colorbar(scatter, ax=ax, label="Protein intensity")
            fig.savefig(intensity_dir_2 / f"Spatial intensity plot - Sample {k} - {marker.replace("*","")}.png", dpi=dpi, facecolor="white")
            plt.close(fig)
    logger.warning(f"Intensities spatial plots generated!")