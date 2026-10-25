import pandas as pd
import numpy as np
from spatial_proteomics.utils import save_anndata_files

import logging
logger = logging.getLogger(__name__)

ABSENT = 0
LOW = 1
MEDIUM = 2
HIGH = 3
MISSING = 255
LABELS = {ABSENT: "Absent", LOW: "Low", MEDIUM: "Medium", HIGH: "High", MISSING: "Missing"}

def _validate_mapping(adata, pm_mapping):
    var_index = pd.Index(adata.var_names)
    result = []
    for pos_idx, (pos_name, intensity_name, pm_short_name) in enumerate(pm_mapping.itertuples(index=False, name=None)):
        if intensity_name not in var_index: continue
        result.append((str(pos_name), str(intensity_name), str(pm_short_name), pos_idx, int(var_index.get_loc(intensity_name))))
        # results contains =  positivity name, intensity name, short name, positivity index, X index
    return result

def _classify_sample(target, pm_mapping):
    """Classify positive cells per marker using *sample-wide positive* P25/P75."""
    mapping = _validate_mapping(target, pm_mapping)
    positivity = np.asarray(target.obsm["Positivity"])
    codes = np.full((target.n_obs, target.n_vars), ABSENT, dtype=np.uint8)
    thresholds = {}
    for pos_name, intensity_name, short_name, p_idx, x_idx in mapping:
        calls = positivity[:, p_idx]
        is_positive = calls == 1
        values = np.asarray(target.X[:, x_idx]).reshape(-1).astype(float)
        valid = is_positive & np.isfinite(values)
        codes[is_positive & ~np.isfinite(values), x_idx] = MISSING
        if valid.any():
            p25, p75 = np.percentile(values[valid], [25, 75])
            codes[valid, x_idx] = MEDIUM
            codes[valid & (values < p25), x_idx] = LOW
            codes[valid & (values > p75), x_idx] = HIGH
        else: p25, p75 = np.nan, np.nan
        thresholds[pos_name] = {"intensity_protein_marker": intensity_name, "pm_short_name": short_name, "p25": float(p25), "p75": float(p75), "n_positive": int(is_positive.sum()), "n_valid": int(valid.sum()),}
    target.obsm["pm_intensity_class"] = codes
    target.uns["pm_intensity_thresholds"] = thresholds
    target.uns["pm_intensity_labels"] = {str(k): v for k, v in LABELS.items()}
    return target

def _summarize_sample(adata, sample_id, ct_short_names):
    """Return (marker summary, observed profiles, phenotype overview) DataFrames.

    Each phenotype uses the positive-marker pattern saved by SpPrAn at
    ``adata.uns[f'sub{short_name}_patterns']``. Only observed profiles are
    enumerated, not all possible combinations. Percentages are percentages
    *within the phenotype*, not percentages of the entire sample.
    """
    codes = np.asarray(adata.obsm["pm_intensity_class"])
    threshold_info = adata.uns["pm_intensity_thresholds"]

    marker_rows, profile_rows, overview_rows = [], [], []
    for short_name in ct_short_names:
        host_col = f"sub{short_name} phenotype"
        patterns_key = f"sub{short_name}_patterns"
        assignments = adata.obs[host_col]
        for phenotype, pattern in adata.uns[patterns_key].items():
            selected = (assignments == phenotype).fillna(False).to_numpy(dtype=bool)
            n_cells = int(selected.sum())
            if n_cells == 0: continue
            positive = [m for m, present in pattern.items() if int(present) == 1]
            columns = []
            for marker in positive:
                intensity_name = threshold_info[marker]["intensity_protein_marker"]
                # pm_short_name = threshold_info[marker]["pm_short_name"]
                x_idx = adata.var_names.get_loc(intensity_name)
                columns.append(codes[selected, x_idx])
                counts = np.bincount(codes[selected, x_idx], minlength=256)
                marker_rows.append({
                    "Sample": sample_id, "Cell type": short_name,
                    "Phenotype": str(phenotype), "Protein marker": marker,
                    "Total cells": n_cells,
                    "Total low intensity cells": int(counts[LOW]), "Total medium intensity cells": int(counts[MEDIUM]),
                    "Total high intensity cells": int(counts[HIGH]),
                    "Low intensity cells proportion (%)": round(counts[LOW] / n_cells * 100,4),
                    "Medium intensity cells proportion (%)": round(counts[MEDIUM] / n_cells * 100,4),
                    "High intensity cells proportion (%)": round(counts[HIGH] / n_cells * 100,4),
                })
            if columns:
                observed = pd.DataFrame(np.column_stack(columns), columns=positive)
                counts_by_profile = observed.value_counts(sort=True, dropna=False)
                for profile_tuple, count in counts_by_profile.items():
                    if len(positive) == 1 and not isinstance(profile_tuple, tuple): profile_tuple = (profile_tuple,)
                    profile_rows.append({
                        "Sample": sample_id, "Cell type": short_name,
                        "Phenotype": str(phenotype),
                        "Intensity profile": " | ".join(f"{threshold_info[marker]["pm_short_name"]}:{LABELS[int(code)]}" for marker, code in zip(positive, profile_tuple)),
                        "Total cells in the profile": int(count), "Percentage (%)": round(int(count) / n_cells * 100,4),
                    })
                n_profiles = len(counts_by_profile)
                top_pct = round(int(counts_by_profile.iloc[0]) / n_cells * 100,4)
            else:
                n_profiles, top_pct = round(0,4), np.nan
            overview_rows.append({
                "Sample": sample_id, "Cell type": short_name,
                "Phenotype": str(phenotype), "Total cells": n_cells,
                "Total of positive protein markers in the phenotype": len(positive),
                "Total of observed profiles": n_profiles,
                "Most common profile proportion (%)": top_pct,
            })
    return (pd.DataFrame(marker_rows), pd.DataFrame(profile_rows), pd.DataFrame(overview_rows))

def phenotype_advanced_cell_classification(adata_dicts, ct_short_names, epm_list, epm_short_names, output_dir, save_anndata, overwrite_existing_files):
    """Classify and summarize samples; export three aggregate CSVs.

    Returns (classified_adatas, marker_df, profile_df, overview_df).
    """
    logger.warning(f"Advanced phenotype classification started...")
    adv_dir = output_dir / "Statistics" / "results" / "Advanced cell classification"
    adv_dir.mkdir(parents=True, exist_ok=True)

    summary_parts  = []
    profile_parts  = []
    overview_parts = []

    pm_mapping = pd.read_csv(output_dir / "results" / "Positivity column names.csv")
    pm_mapping["Short names"] = ""
    for idx in pm_mapping.index:
        if 'DAPI' in pm_mapping.loc[idx,'Positivity column names']: pm_mapping.loc[idx, "Short names"] = "DAPI"
        else: pm_mapping.loc[idx, "Short names"] = epm_short_names[epm_list.index(str(pm_mapping.loc[idx,'Positivity column names']))]

    for k, adata in adata_dicts.items():
        adata = _classify_sample(adata, pm_mapping)
        marker_df, profile_df, overview_df = _summarize_sample(adata, str(k), ct_short_names)
        summary_parts.append(marker_df)
        profile_parts.append(profile_df)
        overview_parts.append(overview_df)
    combine = lambda frames: pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    summary, profiles, overview = map(combine, (summary_parts, profile_parts, overview_parts))

    summary_file = adv_dir / "Protein marker intensity summary.csv"
    profiles_file = adv_dir / "Intensity profiles analysis.csv"
    overview_file = adv_dir / "Phenotype intensity summary.csv"
    
    if summary_file.exists() and not overwrite_existing_files: logger.warning(f"'Protein marker intensity summary.csv' file already exists. Skipping...")
    else: summary.to_csv(summary_file, index=False)
    if profiles_file.exists() and not overwrite_existing_files: logger.warning(f"'Intensity profiles analysis.csv' file already exists. Skipping...")
    else: profiles.to_csv(profiles_file, index=False)
    if overview_file.exists() and not overwrite_existing_files: logger.warning(f"Neighborhood enrichment analysis file already exists. Skipping...")
    else: overview.to_csv(overview_file, index=False)

    adata_path = output_dir / "adata"
    if save_anndata: save_anndata_files(adata_dicts, adata_path)
    logger.warning(f"Advanced phenotype classification finished!")
    return adata_dicts