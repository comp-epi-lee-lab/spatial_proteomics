from spatial_proteomics.unsupervised_mode import _apply_unsupervised_mode
from spatial_proteomics.supervised_mode import _apply_supervised_mode
from spatial_proteomics.anndata import _anndata_generator

import logging
logger = logging.getLogger(__name__)

def _add_phenotype_discovery(adata_dicts, config):
    """
    Convert discovered protein marker patterns into supervised subtype definitions.

    Parameters
    ----------
    adata_dicts : dict of str to anndata.AnnData
        AnnData objects produced by phenotype discovery. Observed phenotype
        protein marker mappings are expected in ``adata.uns``.
    config : dict
        Active SpPrAn configuration. It is updated in place with discovered
        subtype definitions and colors.

    Returns
    -------
    dict
        Updated configuration containing discovered phenotypes in
        ``sub_cell_types``.

    Notes
    -----
    Phenotype definitions are combined across samples for each selected parent
    population. Parent populations with more than ``2**9`` unique phenotypes
    are skipped for the subsequent supervised stage.
    """
    config['spatial_plot'] = {}
    config['spatial_plot']['plot_options'] = 'none'
    config["sub_cell_types"] = {}
    for ct, short_ct in zip(config["cell_types_to_explore"], config["cell_types_to_explore_short_name"]):
        patterns_key = f"sub{short_ct}_patterns"
        combined_patterns = {}
        combined_custom_colors = {}
        for sample_name, adata in adata_dicts.items():
            sample_patterns = adata.uns.get(patterns_key, {})
            combined_patterns.update({host: {marker: int(value) for marker, value in pattern.items()} for host, pattern in sample_patterns.items()})
            combined_custom_colors.update({host:"#ffffff" for host in sample_patterns.keys()})
        logger.warning(f"{ct}: {len(combined_patterns)} unique phenotypes across samples")
        if len(combined_patterns) > 2**9:
            logger.warning(f"WARNING:: Total of unique phenotypes in {ct} is too large for supervised analysis. Skipping supervised analysis for this cell type...")
            continue
        if len(combined_patterns)>0:
            config["sub_cell_types"][ct] = combined_patterns
            config["custom_colors"].update(combined_custom_colors)
    return config


def _apply_full_mode(adata_dicts, config):
    """
    Run phenotype discovery followed by the supervised SpPrAn workflow.

    Parameters
    ----------
    adata_dicts : dict of str to anndata.AnnData
        Initial samples containing major cell-type annotations.
    config : dict
        Validated configuration with ``type_of_analysis`` set to ``"full"``.

    Returns
    -------
    None
        Exploratory and subsequent supervised outputs are written to the
        configured output directory.

    Notes
    -----
    Full mode discovers protein marker-defined phenotypes, converts them into subtype
    rules, regenerates the annotated AnnData objects, and then runs the
    supervised workflow.
    """
    adata_dicts_updated = _apply_unsupervised_mode(adata_dicts, config)
    logger.warning("Unsupervised analysis finished! Starting supervised analysis...")
    config_updated = _add_phenotype_discovery(adata_dicts_updated, config)
    adata_dicts_extended = _anndata_generator(config_updated, config_updated['workspace']['output_dir'])
    _apply_supervised_mode(adata_dicts_extended, config_updated)