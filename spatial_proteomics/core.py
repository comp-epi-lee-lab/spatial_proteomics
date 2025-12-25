from spatial_proteomics.utils import *
import argparse
import sys

def sppran_steps(config):
    """
    Applies common steps in Spatial Proteomics Analaysis (SpPrAn) pipeline.
    
    Parameters
    -----
    config_path : str | Path
        Path to the YAML configuration file.
    """
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

def spatial_proteomics_pipeline(config_path):
    """
    Applies Spatial Proteomics Analaysis (SpPrAn) pipeline.

    Parameters
    -----
    config_path : str | Path
        Path to the YAML configuration file.
    """
    try: config = load_config(config_path)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)
    sppran_steps(config)

def parse_args():
    """
    Parse arguments (like --config) from Command Line (CLI)
    """
    parser = argparse.ArgumentParser(
        description="Run the Spatial Proteomics end-to-end cell typing and quantification analysis pipeline."
    )
    parser.add_argument(
        "--config",
        required=True,
        type=str,
        help="Path to the YAML configuration file."
    )
    return parser.parse_args()

def main():
    """
    Applies Spatial Proteomics Analaysis (SpPrAn) pipeline from Command Line (CLI).
    """
    args = parse_args()
    try: config = load_config(args.config)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)
    sppran_steps(config)

if __name__ == "__main__":
    main()