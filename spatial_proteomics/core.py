from .utils import *

def spatial_proteomics_pipeline():
    config      = load_config()
    adata_dicts = create_or_load_anndata(config)
    plot_spatial(
        adata_dicts              = adata_dicts,
        custom_colors            = config['custom_colors'],
        output_dir               = config['workspace']['output_dir'],
        dpi                      = config['spatial_plot']['dpi'],
        size                     = config['spatial_plot']['scatter_point_size'],
        overwrite_existing_files = config['overwrite_existing_files']
    )
    calculate_cell_proportions(
        adata_dicts              = adata_dicts,
        custom_colors            = config['custom_colors'],
        output_dir               = config['workspace']['output_dir'],
        overwrite_existing_files = config['overwrite_existing_files']
    )