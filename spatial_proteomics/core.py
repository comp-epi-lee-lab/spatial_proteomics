from spatial_proteomics.configuration import load_config
from spatial_proteomics.anndata import create_or_load_anndata
from spatial_proteomics.spatial_plot import plot_spatial
from spatial_proteomics.cell_proportions import calculate_cell_proportions
import argparse
import sys

def sppran_steps(config):
    """
    Run the main SpPrAn analysis steps from a validated configuration.

    Parameters
    ----------
    config : dict
        Validated SpPrAn configuration dictionary. The configuration must
        contain workspace settings, plotting options, cell-type definitions,
        colors, and output-control settings as produced by
        :func:`spatial_proteomics.configuration.load_config`.

    Returns
    -------
    None
        Results are written to the configured output directory.

    Notes
    -----
    This function coordinates AnnData creation or loading, optional spatial
    plotting, and calculation of cell-population counts and proportions.
    """
    adata_dicts = create_or_load_anndata(config)
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
        sub_cell_types           = config['sub_cell_types'] if 'sub_cell_types' in config.keys() else {}
    )

def spatial_proteomics_pipeline(config_path):
    """
    Run the complete SpPrAn pipeline from a YAML configuration file.

    Parameters
    ----------
    config_path : str or pathlib.Path
        Path to the SpPrAn YAML configuration file.

    Returns
    -------
    None
        Analysis outputs are written to the output directory specified in the
        configuration file.

    Raises
    ------
    FileNotFoundError
        If the configuration file or required workspace directories do not
        exist.
    spatial_proteomics.configuration.ConfigValidationError
        If one or more configuration settings are invalid.

    Notes
    -----
    The configuration is validated before any spatial proteomics analysis is
    performed.
    """
    config = load_config(config_path)
    sppran_steps(config)

def parse_args():
    """
    Parse command-line arguments for SpPrAn.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments containing the required ``config`` path.
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
    Run SpPrAn from the command-line interface.

    Notes
    -----
    This is the console entry point used by the ``sppran`` command. Validation
    or analysis errors are written to standard error and terminate the process
    with a non-zero exit code.
    """
    args = parse_args()
    try: 
        config = load_config(args.config)
        sppran_steps(config)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)

def main_gui(config_path):
    """
    Run SpPrAn from the graphical user interface.

    Parameters
    ----------
    config_path : str or pathlib.Path
        Path to the YAML configuration selected in the desktop application.

    Returns
    -------
    None
        Analysis outputs are written to the directory specified by the
        configuration file.
    """
    config = load_config(config_path)
    sppran_steps(config)