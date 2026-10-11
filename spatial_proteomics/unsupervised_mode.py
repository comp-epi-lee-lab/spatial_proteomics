from spatial_proteomics.cell_proportion_matrix import calculate_cell_proportion_matrix
from spatial_proteomics.phenotype_analysis import unsupervised_phenotype_exploration_analysis
from spatial_proteomics.intensity_analysis import calculate_pm_intensity_analysis

def _apply_unsupervised_mode(adata_dicts, config):
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
    if config['type_of_analysis'] == 'full': return adata_dicts_updated
    else: return None