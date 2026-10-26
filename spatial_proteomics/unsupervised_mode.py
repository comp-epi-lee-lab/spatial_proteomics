from spatial_proteomics.cell_proportion_matrix import calculate_cell_proportion_matrix
from spatial_proteomics.phenotype_analysis import unsupervised_phenotype_exploration_analysis
from spatial_proteomics.intensity_analysis import calculate_pm_intensity_analysis, plot_spatial_protein_marker_intensities
from spatial_proteomics.neighborhoods import neighborhood_by_fixed_distance
from spatial_proteomics.advanced_cell_classification import phenotype_advanced_cell_classification

def _apply_unsupervised_mode(adata_dicts, config):
    """
    Run the exploratory protein marker-pattern analysis within selected cell populations.

    Parameters
    ----------
    adata_dicts : dict of str to anndata.AnnData
        Per-sample AnnData objects containing major cell-type annotations.
    config : dict
        Validated unsupervised/full-mode configuration containing selected
        parent populations, exploratory protein markers, optional sample groups, and
        optional intensity markers.

    Returns
    -------
    dict of str to anndata.AnnData or None
        Updated AnnData objects are returned when ``type_of_analysis`` is
        ``"full"`` so discovered phenotypes can continue to the supervised
        workflow. Standard ``"unsupervised"`` mode returns ``None``.

    Notes
    -----
    The core workflow performs exploratory marker-pattern discovery within
    user-defined major cell types. Optional v0.3.1 analyses are executed when
    present/enabled in the configuration:

    * ``intensity_spatial_plot`` generates continuous marker-intensity spatial maps.
    * ``nhood_analysis`` performs fixed-radius neighborhood enrichment,
    interaction-matrix, and co-occurrence analyses.
    * ``advanced_cell_classification`` characterizes discovered phenotypes using
    low/medium/high marker-intensity states.

    These analyses complement rather than replace the original positivity-based
    phenotype discovery.
    """
    output_dir = config['workspace']['output_dir']
    calculate_cell_proportion_matrix(
        adata_dicts              = adata_dicts,
        output_dir               = output_dir,
        custom_colors            = config['custom_colors'],
        focus_ct_name_list       = config['cell_types_to_explore'],
        pie_plot                 = config['want_to_plot_pie_plots'],
        overwrite_existing_files = config['overwrite_existing_files']
    )
    adata_dicts_updated = unsupervised_phenotype_exploration_analysis(
        adata_dicts                            = adata_dicts,
        adata_dir                              = output_dir / "adata",
        exploratory_protein_markers            = config['exploratory_protein_markers'],
        exploratory_protein_markers_short_name = config['exploratory_protein_markers_short_name'],
        cell_types_to_explore                  = config['cell_types_to_explore'],
        cell_types_to_explore_short_name       = config['cell_types_to_explore_short_name'],
        dict_ctrl                              = config['dict_ctrl'],
        dict_trtmnt                            = config['dict_trtmnt'],
        output_dir                             = output_dir,
        overwrite_existing_files               = config['overwrite_existing_files']
    )

    calculate_pm_intensity_analysis(
        adata_dicts               = adata_dicts_updated,
        list_ct_short_name        = config['cell_types_to_explore_short_name'],
        which_pm                  = config['protein_markers_for_intensity'],
        output_dir                = output_dir
    )

    if "intensity_spatial_plot" in config:
        plot_spatial_protein_marker_intensities(
            adata_dicts = adata_dicts_updated,
            output_dir  = output_dir,
            point_size  = config['intensity_spatial_plot']['scatter_point_size'],
            dpi         = config['intensity_spatial_plot']['dpi'],
        )

    if "nhood_analysis" in config:
        adata_dicts_updated = neighborhood_by_fixed_distance(
            adata_dicts              = adata_dicts_updated,
            output_dir               = output_dir,
            save_anndata             = config['locally_save_anndata_files'], 
            overwrite_existing_files = config['overwrite_existing_files'],
            neighborhood_radius_size = config['nhood_analysis']['nhood_radius_size'],
            pixel_conversion_rate    = config['nhood_analysis']['pixel_conversion_rate'],
            max_distance_cooccurrence = config['nhood_analysis']['max_distance_cooccurrence'],
            dpi                      = config['nhood_analysis']['dpi'],
        )

    if config["advanced_cell_classification"] is True:
        adata_dicts_updated = phenotype_advanced_cell_classification(
            adata_dicts              = adata_dicts_updated,
            ct_short_names           = config['cell_types_to_explore_short_name'],
            epm_list                 = config['exploratory_protein_markers'],
            epm_short_names          = config['exploratory_protein_markers_short_name'], 
            output_dir               = output_dir,
            save_anndata             = config['locally_save_anndata_files'], 
            overwrite_existing_files = config['overwrite_existing_files'],
        )

    if config['type_of_analysis'] == 'full': return adata_dicts_updated
    else: return None