from spatial_proteomics.spatial_plot import plot_spatial
from spatial_proteomics.cell_proportions import calculate_cell_proportions

def _apply_supervised_mode(adata_dicts, config):
    """
    Run the supervised SpPrAn analysis workflow.

    Parameters
    ----------
    adata_dicts : dict of str to anndata.AnnData
        Per-sample AnnData objects containing protein marker-defined cell annotations.
    config : dict
        Validated SpPrAn configuration containing supervised spatial-plot
        settings, phenotype colors, output settings, and subtype definitions.

    Returns
    -------
    None
        Spatial plots and cell-population proportion outputs are written to the
        configured output directory.

    Notes
    -----
    Supervised mode assumes that populations of interest are defined a priori
    by the user through protein marker rules.
    """
    if config['spatial_plot']['plot_options'] != 'none':
        plot_spatial(
            adata_dicts              = adata_dicts,
            custom_colors            = config['custom_colors'],
            output_dir               = config['workspace']['output_dir'],
            plot_options             = config['spatial_plot']['plot_options'],
            custom_plot              = config['spatial_plot']['custom_plot'],
            dpi                      = config['spatial_plot']['dpi'],
            size                     = config['spatial_plot']['scatter_point_size'],
            overwrite_existing_files = config['overwrite_existing_files']
        )
    calculate_cell_proportions(
        adata_dicts              = adata_dicts,
        custom_colors            = config['custom_colors'],
        output_dir               = config['workspace']['output_dir'],
        overwrite_existing_files = config['overwrite_existing_files'],
        sub_cell_types           = config['sub_cell_types']
    )