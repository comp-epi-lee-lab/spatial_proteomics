"""
Exploratory marker-defined phenotype analysis for SpPrAn.

This module identifies observed binary protein marker patterns within selected
major cell populations, calculates phenotype abundance, stores phenotype labels
and marker definitions in AnnData objects, and optionally compares phenotype
abundance between two user-defined sample groups.

The workflow enumerates observed marker patterns rather than performing
latent-space clustering.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import ast
import numpy as np
import textwrap
from matplotlib.patches import Patch
from scipy.stats import ttest_ind, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from spatial_proteomics.utils import median_iqr, save_anndata_files
import logging
logger = logging.getLogger(__name__)

def _phenotypes_proportion_matrix(results_dir, subCT_dict, ct_short, overwrite_existing_files=False):
    """
    Export a sample-by-phenotype proportion matrix for one parent population.

    Parameters
    ----------
    results_dir : pathlib.Path
        Base statistics-results directory.
    subCT_dict : dict of str to pandas.DataFrame
        Per-sample phenotype tables containing host names and proportions.
    ct_short : str
        Short name of the parent cell population.
    overwrite_existing_files : bool, default=False
        Whether an existing matrix may be replaced.

    Returns
    -------
    None
        The phenotype-proportion matrix is written as CSV.
    """
    save_filepath = results_dir / "Phenotypes proportion"
    save_filepath.mkdir(parents=True, exist_ok=True)
    save_namefile = save_filepath / f"Phenotypes proportion matrix - sub{ct_short}.csv"
    if save_namefile.exists() and not overwrite_existing_files:
        logger.warning(f"File 'Phenotypes proportion matrix.csv' already exists...")
    else:
        pd.concat([df[['Host name','Proportion']].set_index('Host name').rename(columns={'Proportion':key}).apply(lambda x:x.str.replace('%','').astype(float)) for key,df in subCT_dict.items()], axis=1).to_csv(save_namefile)
        logger.warning(f'Phenotype proportion matrix for {ct_short} saved!')

def _filtering_phenotypes(df_tt,subCT_dict, dictionary, samples_id, subCT_OG_dict, apply_filter=False, zeroes_for_statistics = False):
    """
    Summarize phenotype abundance and detection frequency across a sample group.

    Parameters
    ----------
    df_tt : pandas.DataFrame
        Phenotype table indexed by phenotype host name.
    subCT_dict : dict of str to pandas.DataFrame
        Filtered per-sample phenotype-proportion tables.
    dictionary : dict
        Mapping from sample filenames to the group label.
    samples_id : pandas.DataFrame
        Mapping between SpPrAn sample identifiers and source filenames.
    subCT_OG_dict : dict of str to pandas.DataFrame
        Unfiltered per-sample phenotype tables.
    apply_filter : bool, default=False
        Whether detection-frequency filter columns should be added.
    zeroes_for_statistics : bool, default=False
        Whether non-detected samples contribute zero phenotype abundance.

    Returns
    -------
    pandas.DataFrame
        Group-level phenotype summaries including detected-only and
        zero-inclusive abundance information.

    Notes
    -----
    Conditional summaries use only samples in which the phenotype was detected.
    Zero-inclusive values answer a different biological/statistical question
    and include absence of the phenotype across samples.
    """
    df_tt['Proportion Average'] = 0.0
    df_tt['Standard deviation'] = 0.0
    df_tt['Conditional Proportion Average'] = 0.0
    df_tt['Conditional Standard deviation'] = 0.0
    df_tt['Phenotype'] = "no_name"
    mean_list = []
    std_list = []
    mean_list_cond = []
    std_list_cond = []
    values_list = []
    values_sample_name_list = []
    values_list_cond = []
    values_sample_name_list_cond = []
    frequency_list_cond = []
    how_many_samples = len(dictionary)
    for host in df_tt.index.tolist():
        list_temp = []
        list_temp_cond = []
        list_sample_name_temp = []
        list_sample_name_temp_cond = []
        frequency_cond = []
        for k, df_phenotype in subCT_dict.items():
            if samples_id.loc[k,'Samples filename'] not in dictionary.keys(): continue
            if host in df_phenotype['Host name'].tolist():
                list_temp.append(float(df_phenotype['Proportion'][df_phenotype['Host name']==host].tolist()[0][:-1])/100)
                list_temp_cond.append(float(df_phenotype['Proportion'][df_phenotype['Host name']==host].tolist()[0][:-1])/100)
                list_sample_name_temp.append(k)
                list_sample_name_temp_cond.append(k)
                corresponding_pheno_list = df_phenotype['Markers name'][df_phenotype['Host name']==host].tolist()
                frequency_cond.append('detected')
                df_tt.loc[host, 'Phenotype'] = corresponding_pheno_list[0] if df_tt.loc[host, 'Phenotype'] == "no_name" else df_tt.loc[host, 'Phenotype']
            else:
                if zeroes_for_statistics: # Adding zeroes to no-filters subSets, for statistics
                    list_temp.append(float(0)/100)
                    list_sample_name_temp.append(k)
        mean_list.append(np.mean(list_temp))
        std_list.append(np.std(list_temp))
        mean_list_cond.append(np.mean(list_temp_cond))
        std_list_cond.append(np.std(list_temp_cond))
        values_list.append([val*100 for val in list_temp])
        values_sample_name_list.append([val for val in list_sample_name_temp])
        values_list_cond.append([val*100 for val in list_temp_cond])
        values_sample_name_list_cond.append([val for val in list_sample_name_temp_cond])
        frequency_list_cond.append(f"{len(frequency_cond)} / {how_many_samples}")
    df_tt['Proportion Average'] = mean_list
    df_tt['Standard deviation'] = std_list
    df_tt['Proportions'] = values_list
    df_tt['Conditional Proportions'] = values_list_cond
    df_tt['Conditional Proportion Average'] = mean_list_cond
    df_tt['Conditional Standard deviation'] = std_list_cond
    df_tt['Conditional Sample Frequency'] = frequency_list_cond
    df_tt['Sample name where proportion comes from'] = values_sample_name_list
    if apply_filter:
        how_many_samples = len(dictionary.keys())
        df_tt['Pass 5% filter']  = [True if t/how_many_samples > 0.05 else False for t in df_tt['Counts']]
        df_tt['Pass 10% filter'] = [True if t/how_many_samples > 0.1  else False for t in df_tt['Counts']]
        df_tt['Pass 15% filter'] = [True if t/how_many_samples > 0.15 else False for t in df_tt['Counts']]
        df_tt['Pass 20% filter'] = [True if t/how_many_samples > 0.2  else False for t in df_tt['Counts']]
        df_tt['Pass 25% filter'] = [True if t/how_many_samples > 0.25 else False for t in df_tt['Counts']]
    df_tt = df_tt.sort_values(by="Conditional Proportion Average", ascending=False)
    return df_tt

def _count_phenotypes(subCT_dict, dictionary, results_dir, subCT_OG_dict, apply_filter=False, zeroes_for_statistics=False):
    """
    Count phenotype detections across samples in one experimental group.

    Returns
    -------
    pandas.DataFrame
        Phenotype detection counts and abundance summaries for the group.
    """
    all_phenotypes = []
    samples_id = pd.read_csv(results_dir / 'Samples Id.csv', dtype=str)
    samples_id.set_index('Samples Id', inplace=True)
    for k, df_phenotype in subCT_dict.items():
        if samples_id.loc[k,'Samples filename'] not in dictionary.keys(): continue
        for host in df_phenotype['Host name']:
            all_phenotypes.append(host)
    count_dict = {item:all_phenotypes.count(item) for item in set(all_phenotypes)}
    df_tt = pd.DataFrame(count_dict.values(),index=count_dict.keys(),columns=['Counts'])
    df_tt['Treatment'] = list(dictionary.values())[0]
    df_tt = _filtering_phenotypes(df_tt, subCT_dict, dictionary, samples_id, subCT_OG_dict, apply_filter, zeroes_for_statistics)
    return df_tt

def _common_n_unique_phenotypes(df_phenotypes_ctrl, name_ctrl, df_phenotypes_trtmnt, name_trtmnt, short_ct_name, results_dir, overwrite_existing_files=False, second_filter = None):
    """
    Separate phenotypes into common and group-specific sets.

    Parameters
    ----------
    df_phenotypes_ctrl, df_phenotypes_trtmnt : pandas.DataFrame
        Phenotype summaries for the two comparison groups.
    name_ctrl, name_trtmnt : str
        Display names of the groups.
    short_ct_name : str
        Short name of the parent population.
    results_dir : pathlib.Path
        Statistics-results directory.
    overwrite_existing_files : bool, default=False
        Whether existing outputs may be replaced.
    second_filter : int or None, default=None
        Optional detection-frequency filter percentage.

    Returns
    -------
    pandas.DataFrame
        Phenotypes detected in both groups and their abundance summaries.
    """
    results_dir_to_save = results_dir / f"sub{short_ct_name} split test analysis"
    results_dir_to_save.mkdir(parents=True, exist_ok=True)
    save_namefile_all = results_dir_to_save / f"Common sub{short_ct_name}.csv"
    save_namefile_ctrl = results_dir_to_save / f"Unique sub{short_ct_name} - {name_ctrl}.csv"
    save_namefile_trtmnt = results_dir_to_save / f"Unique sub{short_ct_name} - {name_trtmnt}.csv"

    if second_filter is not None:
        df_phenotypes_ctrl = df_phenotypes_ctrl[df_phenotypes_ctrl[f'Pass {second_filter} filter']==True]
        df_phenotypes_trtmnt = df_phenotypes_trtmnt[df_phenotypes_trtmnt[f'Pass {second_filter} filter']==True]
    together_list = [host for host in df_phenotypes_ctrl.index.tolist() if host in df_phenotypes_trtmnt.index.tolist()]
    list_unique_ctrl = [host for host in df_phenotypes_ctrl.index.tolist() if host not in together_list]
    list_unique_trtmnt = [host for host in df_phenotypes_trtmnt.index.tolist() if host not in together_list]
    
    df_together = df_phenotypes_ctrl.loc[together_list].sort_values(by="Conditional Proportion Average", ascending=False)
    df_together = df_together[["Phenotype","Conditional Proportion Average", 'Conditional Sample Frequency', "Conditional Proportions", "Proportion Average"]].copy()
    df_together = df_together.rename(columns={"Conditional Proportion Average": f"Detected {name_ctrl} %", "Conditional Sample Frequency": f"Detected Sample Frequency in {name_ctrl}", "Conditional Proportions": f"Detected {name_ctrl} Prop Values", "Proportion Average": f"Filtered {name_ctrl} %",})
    df_together[f"Detected {name_trtmnt} %"] = (df_phenotypes_trtmnt.loc[together_list]["Conditional Proportion Average"]*100).round(2)
    df_together[f"Detected {name_trtmnt} Prop Values"] = df_phenotypes_trtmnt.loc[together_list]["Conditional Proportions"]
    df_together[f"Filtered {name_trtmnt} %"] = (df_phenotypes_trtmnt.loc[together_list]["Proportion Average"]*100).round(2)
    df_together[f"Detected {name_ctrl} %"] = (df_together[f"Detected {name_ctrl} %"]*100).round(2)
    df_together[f"Filtered {name_ctrl} %"] = (df_together[f"Filtered {name_ctrl} %"]*100).round(2)
    df_together[f"Detected Difference {name_trtmnt}-{name_ctrl}"] =  df_together[f"Detected {name_trtmnt} %"] - df_together[f"Detected {name_ctrl} %"]
    df_together['Interpretation'] = [f"Increase after {name_trtmnt}" if value>0 else f"Decrease after {name_trtmnt}" for value in df_together[f"Detected Difference {name_trtmnt}-{name_ctrl}"]]
    df_together[f"Detected Sample Frequency in {name_trtmnt}"] = df_phenotypes_trtmnt.loc[together_list]["Conditional Sample Frequency"]
    df_together = df_together[[
        "Phenotype",
        f"Detected {name_ctrl} %",
        f"Detected {name_trtmnt} %",
        f"Detected Difference {name_trtmnt}-{name_ctrl}",
        "Interpretation",
        f"Detected Sample Frequency in {name_ctrl}",
        f"Detected Sample Frequency in {name_trtmnt}",
        f"Detected {name_ctrl} Prop Values",
        f"Detected {name_trtmnt} Prop Values",
        f"Filtered {name_ctrl} %",
        f"Filtered {name_trtmnt} %",
    ]].copy()
    
    df_unique_ctrl = df_phenotypes_ctrl.loc[list_unique_ctrl].sort_values(by="Conditional Proportion Average", ascending=False)
    df_unique_ctrl = df_unique_ctrl[["Phenotype","Conditional Proportion Average", "Proportion Average", "Treatment"]].copy()
    df_unique_ctrl = df_unique_ctrl.rename(columns={"Conditional Proportion Average":"Detected %","Proportion Average":"Filtered %","Treatment":"Condition"})
    df_unique_ctrl["Detected %"] = (df_unique_ctrl["Detected %"]*100).round(2)
    df_unique_ctrl["Filtered %"] = (df_unique_ctrl["Filtered %"]*100).round(2)

    df_unique_trtmnt = df_phenotypes_trtmnt.loc[list_unique_trtmnt].sort_values(by="Proportion Average", ascending=False)
    df_unique_trtmnt = df_unique_trtmnt[["Phenotype","Conditional Proportion Average","Proportion Average","Treatment"]].copy()
    df_unique_trtmnt = df_unique_trtmnt.rename(columns={"Conditional Proportion Average": "Detected %", "Proportion Average": "Filtered %","Treatment":"Condition"})
    df_unique_trtmnt["Detected %"] = (df_unique_trtmnt["Detected %"]*100).round(2)
    df_unique_trtmnt["Filtered %"] = (df_unique_trtmnt["Filtered %"]*100).round(2)
    
    if save_namefile_all.exists() and not overwrite_existing_files: 
        logger.warning(f"File 'Common sub{short_ct_name}.csv' already exists...")
        return pd.read_csv(save_namefile_all, index_col=0)

    df_together.to_csv(save_namefile_all)
    df_together[f"{name_ctrl} Prop Values"] = df_phenotypes_ctrl.loc[together_list].sort_values(by="Conditional Proportion Average", ascending=False)["Proportions"].copy()
    df_together[f"{name_trtmnt} Prop Values"] = df_phenotypes_trtmnt.loc[together_list]["Proportions"]
    df_unique_ctrl.to_csv(save_namefile_ctrl)
    df_unique_trtmnt.to_csv(save_namefile_trtmnt)

    return df_together

def _plotting_table(df_together, name_ctrl, name_trtmnt, pvals, caption_text, short_ct_name, plots_dir, test_name="test", dpi=300, bbox_inches="tight", filter_applied=True, second_filter=None, overwrite_existing_files=False, max_phenotypes_per_plot=12, ctrl_values_column=None, trtmnt_values_column=None):
    """
    Plot paginated two-group phenotype-abundance comparisons.

    Parameters
    ----------
    df_together : pandas.DataFrame
        Common-phenotype comparison table.
    name_ctrl, name_trtmnt : str
        Group labels.
    pvals : array-like
        Nominal p-values shown beside each phenotype.
    caption_text : str
        Figure-caption text.
    short_ct_name : str
        Short name of the parent population.
    plots_dir : pathlib.Path
        Destination directory.
    test_name : str, default="test"
        Statistical-test label used in output filenames.
    dpi : int, default=300
        Figure resolution.
    filter_applied : bool, default=True
        Whether results correspond to filtered phenotype tables.
    second_filter : int or None, default=None
        Optional secondary filter identifier.
    overwrite_existing_files : bool, default=False
        Whether existing figures may be replaced.
    max_phenotypes_per_plot : int, default=12
        Maximum phenotypes displayed per page.
    ctrl_values_column, trtmnt_values_column : str or None
        Columns containing per-sample abundance vectors.

    Returns
    -------
    None
        Plot images and phenotype-ID mappings are written below the
        ``Statistics`` output directory.
    """
    plots_dir.mkdir(parents=True, exist_ok=True)

    if max_phenotypes_per_plot < 1:
        logger.warning("'max_phenotypes_per_plot' must be >= 1. Changing to default value...")
        max_phenotypes_per_plot = 12
    if len(df_together) == 0:
        logger.warning(f"No phenotypes to plot in {short_ct_name}")
        return
    if len(pvals) != len(df_together):
        logger.warning(f"Number of p-values does not match phenotypes in {short_ct_name}. Skipping this plot...")
        return

    # File names
    prefix = f"{second_filter} " if filter_applied and second_filter is not None else ""
    suffix = "" if filter_applied else " - No filter"
    base_name = f"{prefix}Common sub{short_ct_name} proportion differences analysis - Box plot - {test_name}{suffix}"

    # Prepare values
    def to_array(values):
        if isinstance(values, str): values = ast.literal_eval(values)
        if values is None: return np.array([], dtype=float)
        arr = np.asarray(values, dtype=float).ravel()
        return arr[np.isfinite(arr)]

    records = []
    for i, (phenotype, row) in enumerate(df_together.iterrows()):
        ctrl_raw = row[ctrl_values_column] if ctrl_values_column is not None else row.iloc[-4]
        trtmnt_raw = row[trtmnt_values_column] if trtmnt_values_column is not None else row.iloc[-3]
        records.append({"id": f"P{i+1:03d}", "phenotype": str(phenotype), "ctrl": to_array(ctrl_raw), "trtmnt": to_array(trtmnt_raw), "pval": float(pvals[i])})

    # Save phenotype ID mapping
    mapping = pd.DataFrame({"Plot ID": [r["id"] for r in records], "Phenotype": [r["phenotype"] for r in records]})
    mapping_file = plots_dir / f"{base_name} - Phenotype mapping.csv"
    if overwrite_existing_files or not mapping_file.exists(): mapping.to_csv(mapping_file, index=False)

    # Prepare figure pages
    n_pages = int(np.ceil(len(records) / max_phenotypes_per_plot))
    colors = {name_ctrl: "#D9D9D9", name_trtmnt: "#A6C9EC"}
    rng = np.random.default_rng(42)

    def format_pvalue(p):
        if not np.isfinite(p): return "N/A"
        if p < 0.001: return "<0.001"
        return f"{p:.3f}"

    # Generate plots
    for page in range(n_pages):
        start = page * max_phenotypes_per_plot
        end = min(start + max_phenotypes_per_plot, len(records))
        subset = records[start:end]
        n_rows = len(subset)
        
        page_suffix = f" - Page {page+1}" if n_pages > 1 else ""
        save_namefile = plots_dir / f"{base_name}{page_suffix}.png"
        if save_namefile.exists() and not overwrite_existing_files: continue

        # Smaller figures for fewer phenotypes
        fig_height = max(3.5, 0.65 * n_rows + 2.1)
        fig, ax = plt.subplots(figsize=(10, fig_height))
        fig.subplots_adjust(left=0.13, right=0.83, top=0.88, bottom=0.23)
        all_values = []

        # Plot each phenotype
        for j, record in enumerate(subset):
            groups = [(record["ctrl"], j - 0.17, name_ctrl), (record["trtmnt"], j + 0.17, name_trtmnt)]
            for values, position, group in groups:
                if len(values) == 0: continue
                all_values.extend(values.tolist())
                
                if len(values) >= 2: bp = ax.boxplot([values], positions=[position], widths=0.24, vert=False, patch_artist=True, showfliers=False, manage_ticks=False, whis=1.5, medianprops={"color": "#222222", "linewidth": 1.7}, boxprops={"facecolor": colors[group], "edgecolor": "#666666", "linewidth": 0.9}, whiskerprops={"color": "#666666", "linewidth": 1.0}, capprops={"color": "#666666", "linewidth": 1.0})

                # Overlay all biological samples
                jitter = rng.uniform(-0.035, 0.035, len(values))
                ax.scatter(values, position + jitter, s=28, color=colors[group], edgecolor="#555555", linewidth=0.5, alpha=0.85, zorder=5)

            # p-value labels
            ax.text(1.025, j, format_pvalue(record["pval"]), transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=9, clip_on=False)

        # Axis formatting
        ax.set_yticks(range(n_rows))
        ax.set_yticklabels([f"{r['id']} (n={len(r['ctrl'])}/{len(r['trtmnt'])})" for r in subset], fontsize=9)
        ax.set_ylim(n_rows - 0.5, -0.5)
        ax.set_ylabel("Phenotype", fontsize=11)
        ax.set_xlabel("Phenotype proportion (%)", fontsize=11)
        max_value = max(all_values, default=0)
        ax.set_xlim(0, max(1.0, max_value * 1.08))
        ax.grid(axis="x", linestyle="--", alpha=0.25)
        ax.set_axisbelow(True)

        # Legend
        legend_handles = [Patch(facecolor=colors[name_ctrl], edgecolor="#666666", label=name_ctrl), Patch(facecolor=colors[name_trtmnt], edgecolor="#666666", label=name_trtmnt)]
        ax.legend(handles=legend_handles, loc="lower right", fontsize=9, frameon=True)

        # Titles and annotations
        title = f"Conditional phenotype abundance: {name_ctrl} vs {name_trtmnt} in {short_ct_name}"
        if n_pages > 1: title += f" (Page {page+1}/{n_pages})"
        ax.set_title(title, fontsize=13, pad=20)
        ax.text(1.025, 1.02, "Nominal p", transform=ax.transAxes, fontsize=9, fontweight="bold", ha="left")

        # Caption
        wrapped_caption = textwrap.fill(caption_text.replace("\n", " "), width=100)
        fig.text(0.13, 0.045, wrapped_caption, ha="left", va="bottom", fontsize=8)

        # Save
        fig.savefig(save_namefile, dpi=dpi, bbox_inches=bbox_inches)
        plt.close(fig)
        logger.warning(f"Saved boxplot for {short_ct_name}.")

def _common_phenotypes_differences_analysis(df_together, name_ctrl, name_trtmnt, short_ct_name, output_dir, filter_applied=True, second_filter = None, overwrite_existing_files=False):
    """
    Compare common phenotype abundances between two sample groups.

    Parameters
    ----------
    df_together : pandas.DataFrame
        Common-phenotype table containing per-sample abundance vectors.
    name_ctrl, name_trtmnt : str
        Names of the two comparison groups.
    short_ct_name : str
        Short name of the parent population.
    output_dir : pathlib.Path
        SpPrAn output directory.
    filter_applied : bool, default=True
        Whether the analysis uses filtered phenotype tables.
    second_filter : int or None, default=None
        Optional secondary detection-frequency filter.
    overwrite_existing_files : bool, default=False
        Whether existing outputs may be replaced.

    Returns
    -------
    None
        Statistical summaries and comparison plots are written below the
        ``Statistics`` output directory.

    Notes
    -----
    The function performs two-sided Mann–Whitney U and Welch independent-sample
    t-tests. Benjamini–Hochberg FDR-adjusted p-values are also calculated.

    Statistical interpretation should consider biological sample size,
    phenotype detection frequency, effect magnitude, and multiplicity.
    """
    results_dir = output_dir / "Statistics" / "results" / f"sub{short_ct_name} split test analysis"
    plots_dir = output_dir / "Statistics" / "plots" / "phenotypes differences analysis"
    results_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    save_namefile = results_dir / f"Common sub{short_ct_name} proportion differences analysis.csv"
    if save_namefile.exists() and not overwrite_existing_files: 
        logger.warning(f"'proportion differences analysis' csv file already exists...")
        results_table = pd.read_csv(save_namefile, index_col=0)
        pvals   = list(results_table['Mann–Whitney p-value'])
        pvals_m = list(results_table['Raw p-value'])
    else:
        cell_types = df_together.index  # rows are cell types
        pvals = []
        pvals_m = []
        summary = []
        
        for cell in cell_types:
            ctrl_vals = df_together.loc[cell,f"Detected {name_ctrl} Prop Values"]
            treat_vals = df_together.loc[cell,f"Detected {name_trtmnt} Prop Values"]
            
            # mean ± SD
            ctrl_mean, ctrl_sd = np.mean(ctrl_vals), np.std(ctrl_vals, ddof=1)
            treat_mean, treat_sd = np.mean(treat_vals), np.std(treat_vals, ddof=1)
            
            ctrl_med, ctrl_q1, ctrl_q3 = median_iqr(ctrl_vals)
            treat_med, treat_q1, treat_q3 = median_iqr(treat_vals)
        
            # Welch t-test
            sta_m, pval_m = ttest_ind(ctrl_vals, treat_vals, equal_var=False)
            stat, pval = mannwhitneyu(ctrl_vals, treat_vals, alternative="two-sided")
            pvals.append(pval)
            pvals_m.append(pval_m)
            
            summary.append({
                "Cell type": cell,
                "Control (mean ± SD)": f"{ctrl_mean:.2f} ± {ctrl_sd:.2f}",
                "Treatment (mean ± SD)": f"{treat_mean:.2f} ± {treat_sd:.2f}",
                "Control (median [IQR])": f"{ctrl_med:.2f} [{ctrl_q1:.2f}–{ctrl_q3:.2f}]",
                "Treatment (median [IQR])": f"{treat_med:.2f} [{treat_q1:.2f}–{treat_q3:.2f}]",
                "Mann–Whitney p-value": pval,
                "Raw p-value": pval_m
            })
        
        # FDR correction
        if len(pvals)>0: _, fdr_pvals, _, _ = multipletests(pvals, method="fdr_bh")
        else: fdr_pvals = [np.nan]*len(summary)
        if len(pvals_m)>0: _, fdr_pvals_m, _, _ = multipletests(pvals_m, method="fdr_bh")
        else: fdr_pvals_m = [np.nan]*len(summary)
        
        # Add corrected p-values
        for i, row in enumerate(summary):
            row["FDR-corrected M-W p-value"] = fdr_pvals[i]
            row["FDR-corrected RAW p-value"] = fdr_pvals_m[i]
        
        results_table = pd.DataFrame(summary)
        results_table.to_csv(save_namefile)

    _plotting_table(df_together, name_ctrl, name_trtmnt, pvals,   caption_text = "", short_ct_name=short_ct_name, plots_dir=plots_dir, test_name = "Mann–Whitney", filter_applied = filter_applied, second_filter = second_filter, overwrite_existing_files=overwrite_existing_files, ctrl_values_column=f"Detected {name_ctrl} Prop Values", trtmnt_values_column=f"Detected {name_trtmnt} Prop Values") # t test
    _plotting_table(df_together, name_ctrl, name_trtmnt, pvals_m, caption_text = "", short_ct_name=short_ct_name, plots_dir=plots_dir, test_name = "t-test",       filter_applied = filter_applied, second_filter = second_filter, overwrite_existing_files=overwrite_existing_files, ctrl_values_column=f"Detected {name_ctrl} Prop Values", trtmnt_values_column=f"Detected {name_trtmnt} Prop Values") # Mann–Whitney


def _naming_phenotypes(CST_markers, initials_CST, short_ct_name, observed_numbers):
    """
    Create names and marker patterns for observed binary phenotypes.

    Parameters
    ----------
    CST_markers : list of str
        Protein marker columns represented by the binary phenotype code.
    initials_CST : list of str
        Short protein marker names used in readable phenotype labels.
    short_ct_name : str
        Short name of the parent population.
    observed_numbers : array-like of int
        Integer-encoded protein marker patterns observed in the data.

    Returns
    -------
    pandas.DataFrame
        Integer codes, binary patterns, protein marker labels, and generated host names.
    """
    n_markers = len(CST_markers)
    df_bin = pd.DataFrame({"number": sorted(set(observed_numbers))})
    df_bin["binary"] = [format(n, f"0{n_markers}b") for n in df_bin["number"]]
    df_bin["list_pos"] = [list(binary) for binary in df_bin["binary"]]
    df_bin["Markers name"] = [",".join( initials_CST[i] for i, bit in enumerate(binary) if bit == "1") for binary in df_bin["binary"]]
    df_bin['host_name'] = ["".join(["Host_",str(sum([int(s) for s in pos])),"_"]+[initials_CST[n-1] for n in [int(pos[numb_pos])*(numb_pos+1) for numb_pos in range(n_markers)] if n!=0]+[f"_sub{short_ct_name}"]) for pos in df_bin['list_pos']]
    return df_bin

def _calculate_unsupervised_subCT_proportions(adata_dicts, adata_dir, pheno_dict, phenotype_codes_dict, exploratory_protein_markers, output_dir, short_ct_name, ct_name="cells", use_filter=True, overwrite_existing_files=False, min_proportion=0.005):
    """
    Annotate discovered phenotypes and calculate their within-parent abundance.

    Parameters
    ----------
    adata_dicts : dict of str to anndata.AnnData
        Per-sample annotated data.
    adata_dir : pathlib.Path
        Directory used when updated AnnData files are saved.
    pheno_dict : dict of str to pandas.DataFrame
        Per-sample phenotype-code lookup tables.
    phenotype_codes_dict : dict of str to numpy.ndarray
        Integer phenotype code for every cell in each sample.
    exploratory_protein_markers : list of str
        Protein markers represented by the discovered patterns.
    output_dir : pathlib.Path
        SpPrAn output directory.
    short_ct_name : str
        Short name of the parent population.
    ct_name : str, default="cells"
        Parent cell type within which proportions are calculated.
    use_filter : bool, default=True
        Whether low-abundance phenotypes are filtered from exported tables.
    overwrite_existing_files : bool, default=False
        Whether existing phenotype files may be replaced.
    min_proportion : float, default=0.005
        Minimum within-parent proportion retained when filtering.

    Returns
    -------
    subCT_dict : dict of str to pandas.DataFrame
        Filtered per-sample phenotype tables.
    adata_dicts : dict of str to anndata.AnnData
        Updated AnnData objects containing phenotype labels/patterns.
    subCT_OG_dict : dict of str to pandas.DataFrame
        Unfiltered per-sample phenotype tables.

    Notes
    -----
    Proportions are relative to the selected parent cell type, not all cells in
    the tissue sample.
    """
    statistics_results_dir = output_dir / "Statistics" / "results"
    save_filepath = ( statistics_results_dir / f"sub{short_ct_name} cell proportions")
    save_filepath.mkdir(parents=True, exist_ok=True)

    subCT_dict = {}
    subCT_OG_dict = {}

    save_anndata_check = False

    ct_obs_column = f"Only {ct_name}"
    host_column = f"sub{short_ct_name} phenotype"
    patterns_key = f"sub{short_ct_name}_patterns"

    for k, adata in adata_dicts.items():
        df_bin = pheno_dict[k]
        phenotype_codes = np.asarray(phenotype_codes_dict[k])
        
        # One phenotype code per cell
        if len(phenotype_codes) != adata.n_obs:
            logger.warning(f"Sample {k}: phenotype codes and AnnData have different numbers of cells. Skipping this sample...")
            continue
        # Map integer phenotype codes to host names, then to binary marker patterns
        code_to_host = dict(zip( df_bin["number"].astype(np.uint64), df_bin["host_name"]))
        all_host_patterns = dict(zip( df_bin["host_name"], df_bin["list_pos"]))

        if len(exploratory_protein_markers) != len( df_bin["list_pos"].iloc[0]) if not df_bin.empty else False:
            logger.warning(f"Sample {k}: protein marker names and phenotype patterns have different lengths. Skipping this sample...")
            continue
        # Identify cells belonging to the selected cell type
        if ct_obs_column not in adata.obs.columns:
            logger.warning(f"Column '{ct_obs_column}' not found in sample {k}. Skipping this sample...")
            continue

        ct_mask = (adata.obs[ct_obs_column].to_numpy() == ct_name)
        length_CT = int(ct_mask.sum())

        if host_column not in adata.obs.columns: adata.obs[host_column] = pd.Series( pd.NA, index=adata.obs_names, dtype="string")
        else: adata.obs[host_column] = adata.obs[host_column].astype("string")
        adata.obs.loc[ct_mask, host_column] = pd.NA

        if length_CT == 0:
            phenotype_counts = np.zeros( len(df_bin), dtype=np.int64)
            adata.uns[patterns_key] = {}
        else:
            save_anndata_check = True
            ct_codes = phenotype_codes[ct_mask]
            unique_codes, counts = np.unique( ct_codes, return_counts=True)
            count_lookup = dict(zip( unique_codes.tolist(), counts.tolist()))
            phenotype_counts = np.array( [count_lookup.get(int(code), 0) for code in df_bin["number"]], dtype=np.int64)
            cell_host_names = np.array([ code_to_host[int(code)] for code in ct_codes], dtype=object)
            adata.obs.loc[ ct_mask, host_column] = cell_host_names
            observed_hosts = pd.unique(cell_host_names)
            adata.uns[patterns_key] = {host: {marker: value for marker, value in zip( exploratory_protein_markers, all_host_patterns[host])} for host in observed_hosts}

        df_CST = pd.DataFrame({"Host name": df_bin["host_name"].to_numpy(),"Total cells": phenotype_counts})
        if length_CT != 0: proportions_numeric = (phenotype_counts / length_CT)
        else: proportions_numeric = np.zeros( len(phenotype_counts), dtype=float)
        df_CST["Proportion"] = [ f"{p:.4%}" for p in proportions_numeric]
        df_temp = pd.concat([ df_CST.reset_index(drop=True), df_bin["Markers name"].reset_index(drop=True)], axis=1)

        # Apply the original minimum-proportion filter (also known as "first filter")
        if use_filter:
            keep_mask = ((phenotype_counts != 0)& (proportions_numeric >= min_proportion))
            df_to_return = df_temp.loc[keep_mask].copy()
        else:
            df_to_return = df_temp

        title_name = f"Sample {k} ({adata.n_obs} cells)"
        save_namefile_subs = save_filepath / f"sub{short_ct_name} cell proportions - {title_name}.csv"
        if (overwrite_existing_files or not save_namefile_subs.exists()): df_to_return.to_csv(save_namefile_subs)
        subCT_dict[k] = df_to_return
        subCT_OG_dict[k] = df_temp

    if save_anndata_check: save_anndata_files(adata_dicts, adata_dir)

    return subCT_dict, adata_dicts, subCT_OG_dict

def unsupervised_phenotype_exploration_analysis(adata_dicts, adata_dir, exploratory_protein_markers, exploratory_protein_markers_short_name, cell_types_to_explore, cell_types_to_explore_short_name, dict_ctrl, dict_trtmnt, output_dir,overwrite_existing_files):
    """
    Discover and quantify observed protein marker-defined phenotypes.

    Parameters
    ----------
    adata_dicts : dict of str to anndata.AnnData
        Samples containing major cell-type annotations and binary positivity
        data.
    adata_dir : pathlib.Path
        Directory used to save updated AnnData objects.
    exploratory_protein_markers : list of str
        Positivity protein markers selected for phenotype discovery.
    exploratory_protein_markers_short_name : list of str
        Short names corresponding to the exploratory protein markers.
    cell_types_to_explore : list of str
        Major cell populations in which phenotypes are quantified.
    cell_types_to_explore_short_name : list of str
        Short names corresponding to the selected parent populations.
    dict_ctrl, dict_trtmnt : dict or None
        Optional filename-to-group mappings. Both are required for group
        comparison.
    output_dir : pathlib.Path
        SpPrAn output directory.
    overwrite_existing_files : bool
        Whether existing outputs may be replaced.

    Returns
    -------
    dict of str to anndata.AnnData
        Updated AnnData objects containing discovered phenotype labels and
        protein marker-pattern mappings.

    Notes
    -----
    Phenotypes correspond to observed binary positivity combinations rather
    than clusters. When two groups are supplied, common phenotype proportions
    are statistically compared and group-specific phenotypes are reported.
    """
    logger.warning("Started unsupervised_phenotype_exploration_analysis")
    adata_dicts_updated = adata_dicts
    statistics_results_dir = output_dir / "Statistics" / "results" 
    results_og_dir = output_dir / "results" 
    positivity_column_names = pd.read_csv(output_dir / "results" / "Positivity column names.csv")
    positivity_column_names = positivity_column_names[~positivity_column_names["Positivity column names"].str.contains("Positivity - DAPI", na=False)]
    epm_list = positivity_column_names["Name appearing in anndata file"].tolist()

    phenotype_codes_dict = {}
    observed_numbers_dict = {}
    for k, adata in adata_dicts.items():
        phenotype_patterns = pd.DataFrame( adata.obsm["Positivity"], columns=adata.var_names, index=adata.obs_names)[epm_list].to_numpy(dtype=np.uint8)
        n_markers = phenotype_patterns.shape[1]
        if not 1 <= n_markers <= 64: logger.warning(f"Sample {k}: uint64 encoding requires between 1 and 64 markers.")
        if not np.isin(phenotype_patterns, [0, 1]).all(): logger.warning(f"Sample {k}: positivity matrix must contain only binary values.")
        weights = np.uint64(1) << np.arange(n_markers - 1,-1,-1,dtype=np.uint64)
        phenotype_codes = phenotype_patterns.astype(np.uint64) @ weights
        phenotype_codes_dict[k] = phenotype_codes
        observed_numbers_dict[k] = np.unique(phenotype_codes)

    for cell_type_to_explore, cell_type_to_explore_short_name in zip(cell_types_to_explore,cell_types_to_explore_short_name):
        pheno_dict = {}
        for k, adata in adata_dicts_updated.items():
            observed_numbers = observed_numbers_dict[k]
            df_bin = _naming_phenotypes( epm_list, exploratory_protein_markers_short_name, cell_type_to_explore_short_name, observed_numbers)
            pheno_dict[k] = df_bin

        subCT_dict, adata_dicts_updated, subCT_OG_dict = _calculate_unsupervised_subCT_proportions( adata_dicts=adata_dicts_updated, adata_dir=adata_dir, pheno_dict=pheno_dict, phenotype_codes_dict=phenotype_codes_dict, exploratory_protein_markers=exploratory_protein_markers, output_dir=output_dir, short_ct_name=cell_type_to_explore_short_name, ct_name=cell_type_to_explore, use_filter=True, overwrite_existing_files=overwrite_existing_files, min_proportion=0.005)
        _phenotypes_proportion_matrix(results_dir=statistics_results_dir, subCT_dict = subCT_dict, ct_short = cell_type_to_explore_short_name, overwrite_existing_files=overwrite_existing_files)

        if dict_ctrl is not None or dict_trtmnt is not None:
            df_phenotypes_ctrl = _count_phenotypes(subCT_dict, dict_ctrl, results_og_dir, subCT_OG_dict, zeroes_for_statistics=True)
            df_phenotypes_trtmnt = _count_phenotypes(subCT_dict, dict_trtmnt, results_og_dir, subCT_OG_dict, zeroes_for_statistics=True)
            df_together = _common_n_unique_phenotypes(df_phenotypes_ctrl = df_phenotypes_ctrl, name_ctrl = list(dict_ctrl.values())[0], df_phenotypes_trtmnt = df_phenotypes_trtmnt, name_trtmnt = list(dict_trtmnt.values())[0], short_ct_name=cell_type_to_explore_short_name, results_dir=statistics_results_dir, overwrite_existing_files=overwrite_existing_files)
            _common_phenotypes_differences_analysis(df_together = df_together, name_ctrl = list(dict_ctrl.values())[0], name_trtmnt = list(dict_trtmnt.values())[0], short_ct_name=cell_type_to_explore_short_name, output_dir=output_dir, overwrite_existing_files=overwrite_existing_files)

        logger.warning(f"Phenotype analysis for {cell_type_to_explore} completed!")
    logger.warning("All phenotype analyses completed!")

    return adata_dicts_updated