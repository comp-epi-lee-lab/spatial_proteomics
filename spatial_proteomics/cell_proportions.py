# Loading needed libraries

from pathlib import Path
import numpy as np
np.float_ = np.float64
import pandas as pd
import logging
from spatial_proteomics.utils import order_sub_cell_types

logger = logging.getLogger(__name__)

def _safe_percent(numerator, denominator):
    if denominator in (None, 0) or pd.isna(denominator): return 0.0
    return float(np.round((numerator / denominator) * 100, 2))

def _build_hierarchical_report_for_sample(adata, sample_name, custom_colors, sub_cell_types={}, include_other_cells=True):

    rows = []
    total_cells = int(adata.n_obs)
    if total_cells == 0:
        return pd.DataFrame()
    sub_cell_types = sub_cell_types or {}
    parents_n_children = {"General" : [k] for k in [x for x in adata.obs["clusters"].dropna().unique() if x != 'Other cells']} | {parent_key:list(child_dict.keys()) for parent_key, child_dict in sub_cell_types.items()}
    ordered_keys = order_sub_cell_types({"General" : {k:0} for k in [x for x in adata.obs["clusters"].dropna().unique() if x != 'Other cells']} | sub_cell_types)

    def add_row(population, parent, count, denominator):
        rows.append(
            {
                "Sample": sample_name,
                "Population": population,
                "Parent population": parent,
                "Count": int(count),
                "% of parent": _safe_percent(count, denominator),
                "% of total": _safe_percent(count, total_cells),
                # "Denominator": int(denominator) if denominator not in (None, 0) else denominator,
            }
        )

    add_row(population="Total cells", parent="", count=total_cells, denominator=total_cells,)
    for parent in ordered_keys:
        if parent == 'General':
            universe_to_count = adata.obs['clusters'] 
            denominator = adata.n_obs
        else:
            subtype_col = f"Subtypes of {parent}"
            universe_to_count = adata.obs[subtype_col][adata.obs[subtype_col] != "Other cells"]
            denominator = 0
            if parent in adata.obs['clusters'].unique(): denominator = (adata.obs['clusters'] == parent).sum()
            else:
                for idx in adata.obs.columns[adata.obs.columns.str.contains("Subtypes of")].unique():
                    if parent in adata.obs[idx].unique():
                        denominator = (adata.obs[idx] == parent).sum()
                        break
        for child in parents_n_children[parent]:
            add_row(population=child, parent=parent, count=int((universe_to_count == child).sum()), denominator=denominator,)
    return pd.DataFrame(rows)

def calculate_cell_proportions(adata_dicts, custom_colors, output_dir, sub_cell_types={}, overwrite_existing_files=False):
    """
    Quantifies cell population from AnnData and saves it to csv files in results directory.
    
    Parameters
    -----
    adata_dicts : Dict[str:AnnData]
       Dictionary containing multiple annotated data matrices. 
    custom_colors : Dict[str:str]
       Dictionary containing HEX colors per cell type. 
    output_dir : Path
        Path to the output directory.
    overwrite_existing_files : Bool (Optional; Default is False)
        Boolean value to decide if plots will be overwrited in plots directory or not.
    """
    output_dir = Path(output_dir)
    results_dir = output_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    for k, adata in adata_dicts.items():
        if adata.n_obs == 0:
            logger.warning(f"No cells were found in Sample {k}. Skipping sample...")
            continue

        if "color_keys" in adata.uns and "color_vals" in adata.uns:
            for key,val in zip(adata.uns['color_keys'], adata.uns['color_vals']):
                if key not in custom_colors.keys(): custom_colors[key] = val

        for folder_name in adata.uns['folder_names']:
            folder_dir = results_dir / folder_name
            folder_dir.mkdir(parents=True, exist_ok=True)

            if folder_name == 'General':
                universe_to_count = adata.obs['clusters'] 
                total_subcells = adata.n_obs
                title_name = f"Sample {k} ({total_subcells} cells)"
            else:
                subtype_col = f"Subtypes of {folder_name}"
                if subtype_col not in adata.obs.columns:
                    logger.warning(f"{subtype_col} column not found in Sample {k}. Skipping count for this cell type...")
                    continue

                universe_to_count = adata.obs[subtype_col][adata.obs[subtype_col] != "Other cells"]
                total_subcells = 0
                if folder_name in adata.obs['clusters'].unique(): 
                    total_subcells = (adata.obs['clusters'] == folder_name).sum()
                else:
                    for idx in adata.obs.columns[adata.obs.columns.str.contains("Subtypes of")].unique():
                        if folder_name in adata.obs[idx].unique():
                            total_subcells = (adata.obs[idx] == folder_name).sum()
                            break
                if total_subcells == 0: 
                    logger.warning(f"No cells were found for {folder_name} in Sample {k}. Skipping count for this cell type...")
                    continue
                title_name = f"Sample {k} - {folder_name} ({total_subcells} cells)"

            save_namefile = folder_dir / f"Cell type proportions - {title_name}.csv"
            if save_namefile.exists() and not overwrite_existing_files: 
                logger.warning(f"File 'Cell type proportions - {title_name}' already exists. Skipping creation of file...")
                continue

            sum_dict = {"Total cells": int(adata.n_obs)}
            prop_dict = {}

            for cell_type in custom_colors.keys():
                if cell_type not in universe_to_count.values: continue
                sum_dict[cell_type]  = int((universe_to_count == cell_type).sum())
                prop_dict[cell_type] = str(np.round(sum_dict[cell_type]/total_subcells*100,2))+"%"
            
            df = pd.DataFrame(data = {cell_type:[sum_dict[cell_type],prop_dict[cell_type]] for cell_type in custom_colors.keys() if cell_type in universe_to_count.values}, index = ["Total cells in cell type", "Percentage"])
            df.to_csv(save_namefile, index=True)

        report_df = _build_hierarchical_report_for_sample(adata=adata, sample_name=k, custom_colors=custom_colors, sub_cell_types=sub_cell_types)
        report_dir = results_dir / "Cell count report"
        report_dir.mkdir(parents=True, exist_ok=True)
        save_namefile2 = report_dir / f"Cell type proportions report - Sample {k}.csv"
        report_df.to_csv(save_namefile2, index=False)

        print(f"All 'Cell type proportions' csv files for Sample {k} have been created!")

