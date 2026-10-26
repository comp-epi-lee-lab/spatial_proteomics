"""
Fixed-radius spatial neighborhood analyses for SpPrAn.

This module constructs spatial-neighbor graphs from single-cell coordinates and
summarizes relationships among major marker-defined cell populations using
Squidpy neighborhood enrichment, interaction matrices, and co-occurrence
analysis.

Neighborhoods are defined geometrically from a user-specified physical/image
radius rather than from PCA or expression-space nearest neighbors.
"""
import matplotlib.pyplot as plt
import squidpy as sq
import numpy as np
import seaborn as sns

from spatial_proteomics.utils import save_anndata_files

import warnings
import logging
logger = logging.getLogger(__name__)

def _highlighted_by_labels(data_to_analyze, labels, target_row, idx, x_axis_fig_size, y_axis_fig_size, v_bound, rotation, bottom, dpi, title, label_title, save_namefile):
    """
    Save a single-row heatmap highlighting one cell population.

    Parameters
    ----------
    data_to_analyze : numpy.ndarray
        Square neighborhood/interaction matrix.
    labels : pandas.Index or sequence of str
        Cell-population labels corresponding to matrix rows and columns.
    target_row : str
        Cell population being highlighted.
    idx : int
        Row index of ``target_row``.
    x_axis_fig_size, y_axis_fig_size : float
        Figure dimensions in inches.
    v_bound : float
        Symmetric color-scale bound.
    rotation : float
        X-axis label rotation in degrees.
    bottom : float
        Bottom figure margin.
    dpi : int or float
        Output resolution.
    title : str
        Figure title.
    label_title : str
        Color-bar label.
    save_namefile : pathlib.Path
        Destination PNG path.

    Returns
    -------
    None
        The highlighted heatmap is written to disk.
    """
    # idx = labels.get_loc(target_row)
    plt.figure(figsize=(x_axis_fig_size, y_axis_fig_size))
    ax = sns.heatmap(data_to_analyze[[idx], :], xticklabels=labels, yticklabels=[target_row], cmap="viridis", vmin=-v_bound, vmax=v_bound, annot=True, fmt=".2f", cbar_kws={'label': label_title, "ticks": [-v_bound, 0, v_bound]})
    plt.title(title)
    ax.collections[0].colorbar.ax.yaxis.set_label_position("left")
    plt.xticks(rotation=rotation, ha='right')
    fig = plt.gcf()
    fig.set_layout_engine(None)
    fig.subplots_adjust(left=0.05,right=0.95, top=0.85, bottom=bottom)
    fig.savefig(save_namefile,dpi=dpi,bbox_inches=None, pil_kwargs={"compress_level": 1},)
    plt.close(fig)

def _highlight_analysis(adata, sample_idx, save_dir, dpi, overwrite_existing_files,analysis_name, analysis_type=None):
    """
    Generate one highlighted heatmap per major cell population.

    Parameters
    ----------
    adata : anndata.AnnData
        Sample containing neighborhood-analysis results in ``adata.uns`` and
        categorical major cell labels in ``adata.obs["clusters"]``.
    sample_idx : str
        Sample identifier.
    save_dir : pathlib.Path
        Destination directory.
    dpi : int or float
        Output resolution.
    overwrite_existing_files : bool
        Whether existing highlighted images may be replaced.
    analysis_name : str
        Human-readable analysis name used in file names.
    analysis_type : {"nhood", "interaction"} or None
        Matrix to visualize.

    Returns
    -------
    None
        Highlight images are written to disk.
    """
    if analysis_type not in ['nhood','interaction']: return
    if analysis_type=='nhood':
        data_to_analyze = adata.uns["clusters_nhood_enrichment"]["zscore"] 
        h_title = "Neighborhood Enrichment"
        h_label_title='z-score'
    else:
        data_to_analyze = adata.uns["clusters_interactions"] 
        h_title = "Interaction Matrix"
        h_label_title='Interaction proportion'
    labels = adata.obs["clusters"].cat.categories
    v_bound = np.nanmax(np.abs(data_to_analyze))
    max_chars = max([len(lbl) for lbl in labels])
    bottom = min(0.25 + 0.008 * max_chars, 0.75)
    x_axis_fig_size = 5 + 5*((len(labels)-1)//3)
    y_axis_fig_size = 4 if max_chars > 30 else 3
    rotation = 45
    for idx,target_row in enumerate(labels):
        save_namefile = save_dir / f"{analysis_name} - Sample {sample_idx} - {target_row} highlight.png"
        if save_namefile.exists() and not overwrite_existing_files: 
            logger.warning(f"Highlighted {analysis_name} file for {target_row} in sample {sample_idx} already exists. Skipping...")
            continue
        _highlighted_by_labels(data_to_analyze, labels, target_row, idx, x_axis_fig_size, y_axis_fig_size, v_bound, rotation, bottom, dpi, save_namefile=save_namefile, title=f"{h_title} Highlight: {target_row}", label_title=h_label_title)

def _neighborhood_enrichment_analysis(adata,sample_idx,save_dir,dpi,overwrite_existing_files):
    """
    Calculate and plot neighborhood enrichment for one sample.

    Parameters
    ----------
    adata : anndata.AnnData
        Sample with a spatial-neighbor graph and categorical cell labels in
        ``adata.obs["clusters"]``.
    sample_idx : str
        Sample identifier.
    save_dir : pathlib.Path
        Plot-output directory.
    dpi : int or float
        Output resolution.
    overwrite_existing_files : bool
        Whether an existing overview plot may be replaced.

    Returns
    -------
    None
        Squidpy neighborhood-enrichment results are stored in ``adata.uns`` and
        overview/highlight PNG files are written to disk.

    Notes
    -----
    Positive neighborhood-enrichment z-scores indicate cell-population pairs
    observed near one another more often than expected under Squidpy's
    permutation model; negative values indicate spatial avoidance/depletion.
    """
    save_namefile = save_dir / f"Neighborhood enrichment analysis - Sample {sample_idx}.png"
    if save_namefile.exists() and not overwrite_existing_files: 
        logger.warning(f"Neighborhood enrichment analysis file already exists. Skipping...")
        return
    else:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="invalid value encountered in divide",
                category=RuntimeWarning,
                module=r"squidpy\.gr\._nhood"
            )
            sq.gr.nhood_enrichment(adata, cluster_key="clusters", seed=538610, show_progress_bar= False)
        sq.pl.nhood_enrichment(adata, cluster_key="clusters", figsize=(7, 5), show=False)
        fig = plt.gcf()
        fig.set_layout_engine(None)
        fig.savefig(save_namefile,dpi=dpi,bbox_inches=None, pil_kwargs={"compress_level": 1},)
        plt.close(fig)

        _highlight_analysis(adata, sample_idx, save_dir, dpi, overwrite_existing_files, analysis_name='Neighborhood enrichment analysis', analysis_type='nhood')

def _interaction_matrix_analysis(adata,sample_idx,save_dir,dpi,overwrite_existing_files):
    """
    Calculate and plot the cell-population interaction matrix for one sample.

    Parameters
    ----------
    adata : anndata.AnnData
        Sample with a spatial-neighbor graph and ``clusters`` annotations.
    sample_idx : str
        Sample identifier.
    save_dir : pathlib.Path
        Plot-output directory.
    dpi : int or float
        Output resolution.
    overwrite_existing_files : bool
        Whether an existing overview plot may be replaced.

    Returns
    -------
    None
        The interaction matrix is stored in ``adata.uns`` and overview/highlight
        plots are written to disk.
    """
    save_namefile = save_dir / f"Interaction matrix - Sample {sample_idx}.png"
    if save_namefile.exists() and not overwrite_existing_files: 
        logger.warning(f"Interaction matrix analysis file already exists. Skipping...")
        return
    else:
        sq.gr.interaction_matrix(adata, cluster_key="clusters")
        sq.pl.interaction_matrix(adata, cluster_key="clusters", figsize=(7, 5), show=False)
        fig = plt.gcf()
        fig.set_layout_engine(None)
        # fig.subplots_adjust(left=0.22,right=0.86,bottom=0.12,top=0.90)
        fig.savefig(save_namefile,dpi=dpi,bbox_inches=None, pil_kwargs={"compress_level": 1},)
        plt.close(fig)

        _highlight_analysis(adata, sample_idx, save_dir, dpi, overwrite_existing_files, analysis_name='Interaction matrix', analysis_type='interaction')

def _co_occurrence_analysis(adata, sample, save_dir, dpi, overwrite_existing_files, max_distance=None):
    """
    Plot distance-dependent co-occurrence for each major cell population.

    Parameters
    ----------
    adata : anndata.AnnData
        Sample with spatial coordinates, a neighbor graph, and categorical
        ``clusters`` annotations.
    sample : str
        Sample identifier.
    save_dir : pathlib.Path
        Destination directory.
    dpi : int or float
        Output resolution.
    overwrite_existing_files : bool
        Whether existing plots may be replaced.
    max_distance : float or None, default=None
        Optional maximum distance displayed on the x-axis.

    Returns
    -------
    None
        One co-occurrence plot per cell population is written to disk.

    Notes
    -----
    Co-occurrence summarizes how the relative occurrence of population pairs
    changes as a function of spatial distance.
    """
    sq.gr.co_occurrence(adata, cluster_key="clusters")
    for cell_type in adata.obs['clusters'].unique().tolist():
        save_namefile = save_dir / f"Co ocurrence analysis - Sample {sample} - {cell_type}.png"
        if save_namefile.exists() and not overwrite_existing_files:
            logger.warning(f"Co ocurrence analysis file for {cell_type} in sample {sample} already exists...")
            continue
        else:
            sq.pl.co_occurrence(adata, cluster_key="clusters", clusters=cell_type, figsize=(8, 5))
            if max_distance is not None:
                for ax in plt.gcf().axes: ax.set_xlim(right=max_distance)
            plt.savefig(save_namefile,dpi=dpi,bbox_inches=None, pil_kwargs={"compress_level": 1},)
            plt.close()

def neighborhood_by_fixed_distance(adata_dicts, output_dir, dpi, save_anndata, overwrite_existing_files, neighborhood_radius_size=50, pixel_conversion_rate = 1, max_distance_cooccurrence = None):
    """
    Run fixed-radius neighborhood analysis for all samples.

    Parameters
    ----------
    adata_dicts : dict of str to anndata.AnnData
        Per-sample annotated single-cell spatial proteomics data.
    output_dir : pathlib.Path
        SpPrAn output directory.
    dpi : int or float
        Resolution of exported neighborhood plots.
    save_anndata : bool
        Whether AnnData objects containing neighborhood results should be saved.
    overwrite_existing_files : bool
        Whether existing plot files may be replaced.
    neighborhood_radius_size : float, default=50
        User-facing neighborhood radius.
    pixel_conversion_rate : float, default=1
        Conversion factor between the user-facing radius unit and the spatial
        coordinate unit stored in AnnData.
    max_distance_cooccurrence : float or None, default=None
        Optional maximum distance shown in co-occurrence plots.

    Returns
    -------
    dict of str to anndata.AnnData
        AnnData objects containing the computed spatial-neighbor graph and
        Squidpy neighborhood-analysis results.

    Notes
    -----
    The radius supplied to Squidpy is calculated as
    ``neighborhood_radius_size / pixel_conversion_rate``.

    This analysis defines neighborhoods from tissue/image coordinates. It does
    not require PCA, UMAP, or expression-space nearest-neighbor construction.
    """
    logger.warning(f"Neighborhood analysis started...")
    radius_conversion = neighborhood_radius_size / pixel_conversion_rate
    save_dir = output_dir / "Statistics" / "plots" / "Neighborhood analysis"
    save_dir.mkdir(parents=True, exist_ok=True)

    for k, adata in adata_dicts.items():
        save_dir_2 = save_dir / f"{k}"
        save_dir_2.mkdir(parents=True, exist_ok=True)
        sq.gr.spatial_neighbors(adata,radius=radius_conversion, coord_type = 'generic')
        _neighborhood_enrichment_analysis(adata,k,save_dir = save_dir_2,dpi=dpi, overwrite_existing_files=overwrite_existing_files)
        _interaction_matrix_analysis(adata,k,save_dir = save_dir_2,dpi=dpi, overwrite_existing_files=overwrite_existing_files)
        _co_occurrence_analysis(adata,k,save_dir = save_dir_2, dpi=dpi, overwrite_existing_files=overwrite_existing_files, max_distance = max_distance_cooccurrence)

    adata_path = output_dir / "adata"
    if save_anndata: save_anndata_files(adata_dicts,adata_dir=adata_path)
    logger.warning(f"Neighborhood analysis finished!")
    return adata_dicts