import pandas as pd

def _read_intensities(adata_dicts, which_pm = 'all'):
    """
    Extract protein marker intensity matrices from AnnData objects.

    Parameters
    ----------
    adata_dicts : dict of str to anndata.AnnData
        Per-sample AnnData objects.
    which_pm : "all" or list of str, default="all"
        AnnData variable names to retain.

    Returns
    -------
    dict of str to pandas.DataFrame
        Per-sample intensity matrices indexed by cell identifier.
    """
    intensities_dict = {}
    for k, adata in adata_dicts.items():
        if which_pm == 'all': intensities_dict[k] = pd.DataFrame(adata.X,columns=adata.var_names,index=adata.obs_names)
        else: intensities_dict[k] = pd.DataFrame(adata.X,columns=adata.var_names,index=adata.obs_names)[which_pm]
    return intensities_dict

def _boxplot_calculation(df,pm_name):
    """
    Calculate Tukey box plot summary values for one protein marker.

    Returns
    -------
    list of float
        Lower whisker, Q1, median, Q3, and upper whisker using the 1.5-IQR rule.
    """
    q1 = df[pm_name].quantile(0.25)
    median = df[pm_name].quantile(0.50)
    q3 = df[pm_name].quantile(0.75)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    inner_data = df[(df[pm_name] >= lower_bound) & (df[pm_name] <= upper_bound)][pm_name]
    lower_whisker = inner_data.min()
    upper_whisker = inner_data.max()
    return [lower_whisker, q1, median, q3, upper_whisker]

def _boxplot_outliers(df,pm_name):
    """
    Return intensity observations outside the 1.5-IQR whisker range.

    Returns
    -------
    pandas.Series
        Sorted outlier intensity values.
    """
    q1 = df[pm_name].quantile(0.25)
    q3 = df[pm_name].quantile(0.75)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    outliers = df[(df[pm_name] < lower_bound) | (df[pm_name] > upper_bound)][pm_name].sort_values()
    return pd.Series(outliers)

def _intensities_data(intensities_dict,adata_dicts,list_ct_short_name,output_dir):
    """
    Export box plot summary data by discovered phenotype.

    Parameters
    ----------
    intensities_dict : dict of str to pandas.DataFrame
        Per-sample protein marker intensity matrices.
    adata_dicts : dict of str to anndata.AnnData
        AnnData objects containing discovered phenotype labels.
    list_ct_short_name : list of str
        Short names of explored parent populations.
    output_dir : pathlib.Path
        Base directory for intensity analysis results.

    Returns
    -------
    None
        Per-sample/per-phenotype summaries and outlier tables are written as CSV.
    """
    intensities_csv_dir_phen = output_dir / "results" / "Intensity analysis" / "By phenotype"
    intensities_csv_dir_phen.mkdir(parents=True, exist_ok=True)

    for k, adata in adata_dicts.items():
        intensity_df = intensities_dict[k]
        pm_list = intensity_df.columns.tolist()

        for short_ct_name in list_ct_short_name:
            host_column = f"sub{short_ct_name} phenotype"
            
            list_host_names = adata.obs[host_column].unique().tolist()
            for phenotype in list_host_names:
                if pd.isna(phenotype):
                    continue
                title_name_1 = f"Box plot data - Sample {k} - {host_column} - {phenotype}"
                save_namefile_1 = intensities_csv_dir_phen / f"{title_name_1}.csv"
                title_name_2 = f"Box plot data - Sample {k} - {host_column} - {phenotype} - outliers"
                save_namefile_2 = intensities_csv_dir_phen / f"{title_name_2}.csv"

                cell_ids = adata.obs_names[adata.obs[host_column] == phenotype].tolist()
                filtered_df = intensity_df[intensity_df.index.isin(cell_ids)]

                pd.DataFrame({'Metric': ['Lower_Whisker', 'Q1', 'Median (Q2)', 'Q3', 'Upper_Whisker']}|{ pm_name: _boxplot_calculation(filtered_df,pm_name) for pm_name in pm_list }).to_csv(save_namefile_1, index=False)
                pd.DataFrame({ pm_name: _boxplot_outliers(filtered_df,pm_name) for pm_name in pm_list }).to_csv(save_namefile_2, index=False)

def calculate_pm_intensity_analysis(adata_dicts, list_ct_short_name, output_dir, which_pm):
    """
    Export protein marker intensity summaries for discovered phenotypes.

    Parameters
    ----------
    adata_dicts : dict of str to anndata.AnnData
        AnnData objects containing discovered phenotype annotations.
    list_ct_short_name : list of str
        Short names of explored parent populations.
    output_dir : pathlib.Path
        SpPrAn output directory.
    which_pm : list of str
        Positivity column names whose corresponding intensity measurements
        should be summarized. An empty list disables the analysis.

    Returns
    -------
    None
        Summary tables are written under the statistics/intensity output tree.

    Notes
    -----
    This is descriptive output for box plot visualization; the function itself
    does not perform inferential tests on marker intensity.
    """
    if len(which_pm)!=0:
        positivity_column_names = pd.read_csv(output_dir / "results" / "Positivity column names.csv")
        epm_list = positivity_column_names[positivity_column_names["Positivity column names"].isin(which_pm)]["Name appearing in anndata file"].tolist()
        intensities_dict = _read_intensities(adata_dicts, epm_list)
        _intensities_data(intensities_dict, adata_dicts, list_ct_short_name, output_dir / "Statistics")