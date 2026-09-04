import pandas as pd
import scanpy as sc
from anndata import AnnData
from spatial_proteomics.utils import plot_general_cell_hierarchy, plot_digraph, order_sub_cell_types
import logging

logger = logging.getLogger(__name__)

def load_anndata_files(output_dir):
    """
    Load previously generated AnnData objects for all known samples.

    Parameters
    ----------
    output_dir : pathlib.Path
        SpPrAn output directory containing the ``results`` and ``adata``
        subdirectories.

    Returns
    -------
    dict of str to anndata.AnnData or None
        Mapping from sample identifier to loaded AnnData object. A missing
        sample file is represented by ``None``. An empty dictionary is
        returned when the sample index file is unavailable.

    Notes
    -----
    Sample identifiers are read from ``results/Samples Id.csv`` and AnnData
    files are expected at ``adata/Sample_<sample>.h5ad``.
    """
    filepath = output_dir / "results" / "Samples Id.csv"
    if filepath.exists():
        sample_names = pd.read_csv(filepath, dtype=str)['Samples Id'].tolist()
        adata_dicts = {}
        for sample in sample_names:
            sample = str(sample)                                                                                             # Ensuring we can deal with numbers as identifiers.
            adata_path = output_dir / "adata" / f"Sample_{sample}.h5ad"
            if adata_path.exists():
                adata_dicts[sample] = sc.read_h5ad(adata_path)
            else:
                logger.warning(f"WARNING: Anndata loading cannot be completed. Sample_{sample}.h5ad doesn't exist or it's missing.") # An None value is saved in the dictionary.
                adata_dicts[sample] = None
        return adata_dicts
    logger.warning(f"WARNING: Anndata loading cannot be completed. File '{filepath}' doesn't exist or it's missing.")                # An empty dictionary is returned.
    return {}

def filenames(input_dir, filetype):
    """
    Find SpPrAn object-level input files.

    Parameters
    ----------
    input_dir : pathlib.Path
        Directory containing Visiopharm object-level exports.
    filetype : {"tsv", "csv"}
        Input-file extension.

    Returns
    -------
    list of pathlib.Path
        Files matching ``*_objects.<filetype>``.

    Raises
    ------
    FileNotFoundError
        If ``input_dir`` does not exist.

    Notes
    -----
    Hidden files and macOS ``._`` resource-fork files are excluded.
    """
    if input_dir.exists():
        entries = list(input_dir.iterdir())
        filenames = [entry for entry in entries if entry.is_file() and entry.match(f"*_objects.{filetype}") and not entry.name.startswith("._") and not entry.name.startswith(".")]  # Exclude hidden files and macOS resource forks
        return filenames
    else:
        raise FileNotFoundError(f"Input directory '{input_dir}' doesn't exists or it's misspelled.")

def cleaned_data(file_names, output_dir, protein_markers, cell_subtype_dict, filetype='tsv'):
    """
    Clean object-level spatial proteomics tables and create AnnData objects.

    Parameters
    ----------
    file_names : list of pathlib.Path
        Object-level input files to process.
    output_dir : pathlib.Path
        Directory where SpPrAn result metadata will be written.
    protein_markers : list of str
        Primary positivity-column names used for cell-type classification.
    cell_subtype_dict : dict
        Hierarchical subtype definitions. Subtype marker columns are retained
        when required for downstream classification.
    filetype : {"tsv", "csv"}, default="tsv"
        Input table format.

    Returns
    -------
    data_dicts : dict of str to pandas.DataFrame
        Cleaned per-sample tables used for phenotype assignment.
    adata_dicts : dict of str to anndata.AnnData
        Per-sample AnnData objects containing cell-level measurements and
        spatial information.

    Notes
    -----
    When the DAPI positivity column is available, rows lacking nuclear DAPI
    positivity are removed. SpPrAn also records sample identifiers and the
    detected positivity-column names in the output directory.

    The resulting AnnData objects represent individual segmented cells, not
    aggregated tissue regions.
    """
    adata_dicts    = {}
    data_dicts     = {}
    filename_dicts = {}
    if len(file_names)!=0:
        for filename in file_names:
            separator = '\t' if filetype=='tsv' else ','
            filename = str(filename)
            data = pd.read_csv(f"{filename}",sep=separator)
            if data.columns.isin(['Positivity - DAPI (MV - NUC)']).any():
                data = data[data['Positivity - DAPI (MV - NUC)']==1]                                                           # filter out cells without nucleus
            else: logger.warning("WARNING: This pipeline uses 'Positivity - DAPI' for filtering out cells without nucleus. However, such column was not found. Continue without filtering...")
            name_of_file = str(filename.split('/')[-1].split('\\')[-1])
            if name_of_file[:11]=='COMET_6x6 (':
                data['Name'] = [f"{s[11]}{s[14]}_{i:05}" for s,i in zip([name_of_file]*len(data),data.index)]   # retain letter and number for identifying samples
            else:
                data['Name'] = [f"{s}_{i:05}" for s,i in zip([name_of_file[:-12]]*len(data),data.index)]        # retain filename for identifying samples
            data.rename(columns={"Name":"cellID"}, inplace=True)
            
            pos_cols = data.columns[(data.columns.str.contains('Positivity'))&(data.columns.str.contains('MV'))]               # identify columns that says if protein marker is present (1) or abscent (0)
            if len(pos_cols) == 0: 
                logger.warning(f"WARNING: No protein marker valid name was found in Sample {name_of_file[-14:-12]}. Skipping this sample from analysis...")
                continue
            pm_extended = list(set(protein_markers + [pm_e for v1_temp in cell_subtype_dict.values() for v2_temp in v1_temp.values() for pm_e in v2_temp.keys()]))
            for pm_test in pm_extended:
                if pm_test not in pos_cols: logger.warning(f"WARNING: Protein marker '{pm_test}' was not found in Sample {name_of_file[-14:-12]}. Skipping this sample from analysis...")
            if not set(pm_extended).issubset(pos_cols): continue

            data = data[~(data[pos_cols].eq(-1).any(axis=1))]                                                                  # filter out those cells that does contain a NaN value (-1) in the previous columns
            columns_nuc = data.columns[                                                                                        # keep columns with protein intensity in the nucleus
                (data.columns.isin(['cellID']))|
                ((data.columns.str.contains("MV - NUC - "))&(~data.columns.str.contains("Type")))
            ]
            columns_cyto = data.columns[                                                                                        # keep columns with protein intensity in the cytoplasm
                (data.columns.isin(['cellID']))|
                ((data.columns.str.contains("MV - CYTO - "))&(~data.columns.str.contains("Type")))
            ]
            columns_cell = data.columns[                                                                                        # keep columns with protein intensity in the cell
                (data.columns.isin(['cellID']))|
                ((data.columns.str.contains("MV Cell - "))&(~data.columns.str.contains("Type")))
            ]
            intensities_columns = [new_col_check
                for column_check in pos_cols
                for val_check,where_check in {'CYTO':columns_cyto,'NUC':columns_nuc,'Cell':columns_cell}.items()
                if val_check in column_check
                for new_col_check in where_check
                if column_check.replace("DAPI", "DAPI,").replace("Positivity - ",",").replace("*","*,").replace(" (MV", ", (MV").split(",")[1] == new_col_check.rsplit(" - ", 1)[-1].strip()
            ]
            if len(intensities_columns)!=len(pos_cols):
                for column_check,new_col_check in zip(pos_cols,intensities_columns):
                    if column_check.replace("DAPI", "DAPI,").replace("Positivity - ",",").replace("*","*,").replace(" (MV", ", (MV").split(",")[1] != new_col_check.rsplit(" - ", 1)[-1].strip():
                        intensities_columns.insert(pos_cols.tolist().index(column_check),column_check)

            columns_pos = data.columns[data.columns.isin(['cellID', 'X-coordinate', 'Y-coordinate'])]                          # keep spatial columns
            data_temp = data[intensities_columns]#.iloc[:,1:]
            adata = AnnData(                                                                                                   # generate AnnData file to include spatial data
                data_temp.set_index((str(x) for x in data_temp.index)),
                obsm={
                    "spatial": data[columns_pos].iloc[:,1:].to_numpy(),
                    "ID_cell":data[['cellID']].to_numpy(),
                    "Positivity":data[pos_cols].to_numpy()
                },
                uns={
                    "spatial":{"unique":{}},
                    "filename_id":{name_of_file: name_of_file[-14:-12]}
                }
            )

            data = data[["cellID"]+pos_cols.tolist()]
            adata_dicts[f"{data.iloc[0,0][-8:-6]}"] = adata
            data_dicts[f"{data.iloc[0,0][-8:-6]}"] = data
            filename_dicts[f"{data.iloc[0,0][-8:-6]}"] = name_of_file
        df = pd.DataFrame({"Samples Id":filename_dicts.keys(), "Samples filename":filename_dicts.values()})
        results_path = output_dir / "results"
        results_path.mkdir(parents=True, exist_ok=True)
        df.to_csv(results_path / "Samples Id.csv", index=False)
        
        df = pd.DataFrame({"Positivity column names": list(data[pos_cols].columns)})
        df.to_csv(results_path/ "Positivity column names.csv", index=False)
    else: 
        data_dicts = {}
        adata_dicts = {}
    return data_dicts, adata_dicts

def create_cell_type_dict(protein_markers, cell_types):
    """
    Convert ordered cell-type rules into marker-to-state dictionaries.

    Parameters
    ----------
    protein_markers : list of str
        Ordered positivity-column names.
    cell_types : dict of str to list
        Cell-type definitions containing ``1`` for positive, ``0`` for
        negative, and ``None`` for unconstrained markers.

    Returns
    -------
    dict of str to dict
        Marker-based phenotype rules keyed by cell-type name.

    Notes
    -----
    Unconstrained markers are omitted from the returned rule dictionary.
    """
    cell_type_dict = {}
    for cell_type, rule in cell_types.items():
        cell_type_dict[cell_type] = {protein_markers[i]:rule[i] for i in range(len(rule)) if rule[i] is not None}
    return cell_type_dict





def get_index_nones(cell_types, protein_markers, sub_protein_markers, sub_cell_types):
    """
    Identify unconstrained marker positions for each phenotype.

    Parameters
    ----------
    cell_types : dict
        Primary marker-based phenotype definitions.
    protein_markers : list of str
        Primary marker names.
    sub_protein_markers : dict
        Subtype marker names grouped by parent population.
    sub_cell_types : dict
        Hierarchical subtype definitions.

    Returns
    -------
    dict of str to list of int
        Marker positions that are not constrained for each phenotype.
    """
    none_dicc = {}
    total_markers = len(protein_markers)
    # none_subindices = [index+total_markers for index in range(len(sub_protein_markers['all']))]
    for cell_type in cell_types.keys():
        none_indices = [index for index, value in enumerate(protein_markers) if value not in cell_types[cell_type].keys()]
        none_dicc[cell_type]=none_indices
    for cell_type in sub_cell_types.keys():
        # none_indices_temp = none_dicc[cell_type]+none_subindices
        for subcell_type in sub_cell_types[cell_type].keys():
            # none_indices = [index for index in none_indices_temp if index not in [total_markers + sub_protein_markers['all'].index(val) for val in sub_cell_types[cell_type][subcell_type].keys()]]
            none_indices = [total_markers + sub_protein_markers['all'].index(val) for val in sub_protein_markers[cell_type] if val not in sub_cell_types[cell_type][subcell_type].keys()]
            none_dicc[subcell_type]=none_indices
    return none_dicc

def fill_n_order_df_by_total_none_values(df, idx_name, order_dicc):
    """
    Expand unconstrained marker states and order phenotype columns.

    Parameters
    ----------
    df : pandas.DataFrame
        Marker-by-phenotype rule table.
    idx_name : str
        Name assigned to the DataFrame index.
    order_dicc : dict
        Phenotype ordering grouped by specificity/dependency.

    Returns
    -------
    expanded_df : pandas.DataFrame
        Rule table in which unconstrained states are represented as both
        possible binary states.
    original_df : pandas.DataFrame
        Ordered copy retaining the original unconstrained values.
    """
    df_og = df.copy()
    for col in df.columns:
        mask = df[col].isna()
        df.loc[mask, col] = pd.Series( [[0, 1]] * mask.sum(), index=df.index[mask], dtype=object)
    df.index.name = idx_name
    df_og.index.name = idx_name
    df = df[list(reversed([val for cols in order_dicc.values() for val in cols if val in df.columns]))]
    df_og = df_og[list(reversed([val for cols in order_dicc.values() for val in cols if val in df_og.columns]))]
    return df, df_og

def sets_and_subsets_cell_types(cell_type_dict, sub_cell_type_dict, protein_markers, output_dir, custom_colors, arrows_orientation="contains"):
    """
    Determine hierarchical subset relationships among marker-defined phenotypes.

    Parameters
    ----------
    cell_type_dict : dict
        Primary marker-based cell-type definitions.
    sub_cell_type_dict : dict
        User-defined hierarchical subtype definitions.
    protein_markers : list of str
        Primary marker names.
    output_dir : pathlib.Path
        Output directory used for hierarchy visualizations and related files.
    custom_colors : dict of str to str
        Hexadecimal color assigned to each phenotype.
    arrows_orientation : {"contains", "isContained"}, default="contains"
        Direction used to display relationships in hierarchy diagrams.

    Returns
    -------
    final_subsets : dict
        Final subset relationships for each hierarchy level.
    subtype_definitions : dict
        Marker definitions associated with resolved phenotype subsets.
    larger_groups : dict
        Parent or less-specific phenotype groups identified during hierarchy
        resolution.
    custom_colors : dict
        Color mapping including any colors generated for intersecting
        phenotypes.

    Notes
    -----
    Cell populations are compared as sets defined by marker positivity rules.
    These subset relationships describe the logical phenotype definitions used
    by SpPrAn and should not automatically be interpreted as developmental
    lineage relationships.
    """
    sub_protein_markers = {cell_type: [] for cell_type in sub_cell_type_dict.keys()}
    all_temp = []
    for cell_types in sub_cell_type_dict.keys():
        sublist_temp = [subprotein_marker for sub_cell_type in sub_cell_type_dict[cell_types].keys() for subprotein_marker in sub_cell_type_dict[cell_types][sub_cell_type].keys()]
        sub_protein_markers[cell_types] = list(set(sublist_temp))
        all_temp.extend(sub_protein_markers[cell_types])
    sub_protein_markers['all'] = list(set(all_temp))
    none_dicc = get_index_nones(cell_type_dict, protein_markers, sub_protein_markers, sub_cell_type_dict)
    
    # Reordering cell types based on the number of None values in their definitions
    max_N = max([len(N_pos) for N_pos in none_dicc.values()])
    dicc_temp = {}
    for i in reversed(range(max_N+1)): dicc_temp[i] = []
    for cell_type, N_pos in none_dicc.items(): dicc_temp[len(N_pos)].append(cell_type)

    # Creating DataFrames for cell types and sub cell types
    ct_df = pd.DataFrame(cell_type_dict, index=protein_markers, dtype=object)
    ct_df, ct_df_og = fill_n_order_df_by_total_none_values(ct_df, 'Protein Marker', dicc_temp)
    sct_dfs = {}
    # sct_dfs_og = {}
    for sct, data in sub_cell_type_dict.items():
        if sct in cell_type_dict.keys():
            full_index = protein_markers + [p for p in sub_protein_markers[sct] if p not in protein_markers]
            df = pd.DataFrame(data, index=full_index, dtype=object)
            for prot in protein_markers: df.loc[prot] = cell_type_dict[sct].get(prot, None)
        else:
            whole_protein_markers = {}
            hold_temp = [sct]
            while len(hold_temp)>0:
                for key in hold_temp:
                    if key in cell_type_dict.keys():
                        whole_protein_markers.update(cell_type_dict[key])
                        hold_temp.remove(key)
                        break
                    for k in sub_cell_type_dict.keys():
                        if k==key: continue
                        if key in sub_cell_type_dict[k].keys():
                            hold_temp.append(k)
                            whole_protein_markers.update(sub_cell_type_dict[k][key])
                    hold_temp.remove(key)
            full_pm = list(set(whole_protein_markers.keys()))
            full_index = full_pm + [p for p in sub_protein_markers[sct] if p not in full_pm]
            df = pd.DataFrame(data, index=full_index, dtype=object)
            for prot in full_pm: df.loc[prot] = whole_protein_markers.get(prot, None)
        df, df_og = fill_n_order_df_by_total_none_values(df, 'Protein Marker', dicc_temp)
        sct_dfs[sct] = df
        # sct_dfs_og[sct] = df_og

    list_final_subsets = {}
    list_sub_cell_types = {}
    list_bigger_cell_types = {}
    plot_general_cell_hierarchy(protein_markers,cell_type_dict,sub_cell_type_dict,custom_colors,save_path = output_dir / "plots" , arrows_orientation=arrows_orientation)
    # raise Exception ("clJust need to see the new plot")
    list_final_subsets["General"], list_sub_cell_types["General"], list_bigger_cell_types["General"], custom_colors = plot_digraph(ct_df, protein_markers, output_dir, "General", dicc_temp, custom_colors) 
    for sct, df in sct_dfs.items():
        list_final_subsets[sct], list_sub_cell_types[sct], list_bigger_cell_types[sct], custom_colors = plot_digraph(df, list(df.index), output_dir, sct, dicc_temp, custom_colors)
    return list_final_subsets, list_sub_cell_types, list_bigger_cell_types, custom_colors

def assign_cell_types_vectorized(data, cell_type_dict, container):
    """
    Assign marker-defined phenotype labels to cells using vectorized comparisons.

    Parameters
    ----------
    data : pandas.DataFrame
        Cell-level table containing the positivity columns referenced by the
        phenotype rules.
    cell_type_dict : dict
        Mapping from phenotype names to marker-state rules.
    container : str
        Name of the parent population being subdivided. Use ``"General"`` for
        primary cell typing.

    Returns
    -------
    pandas.Series
        Phenotype label for every row in ``data``.

    Notes
    -----
    Cells are assigned in dictionary iteration order and only previously
    unassigned cells are evaluated for subsequent phenotypes. Cells not
    matching any rule are labeled ``"Other cells"`` at the general level.
    """
    labels = pd.Series("Other cells", index=data.index)
    unassigned = pd.Series(True, index=data.index)
    for cell_type, rule in cell_type_dict.items():
        mask = unassigned.copy()
        for marker, expected in rule.items():
            if expected is None: continue
            if isinstance(expected, (list, tuple, set)): mask &= data[marker].isin(expected)
            else: mask &= data[marker] == expected
        labels[mask] = cell_type if cell_type!=container else f"Other {cell_type}"
        unassigned &= ~mask
    return labels

def get_all_children(celltype, dicc, visited=None):
    """
    Return all descendants of a cell type from a nested hierarchy.

    Parameters
    -----
    celltype : str
        The cell type for which to find all children.
    dicc : Dict[str:List]
        Dictionary containing cell types and their corresponding children.
    visited : set (Optional)
        Set of already visited cell types to avoid infinite recursion. 

    Returns
    ------
    set
       Set containing all children of the given cell type.
    """
    if visited is None:
        visited = set()
    
    if celltype in visited:
        return set()
    
    visited.add(celltype)
    node = dicc.get(celltype, [])
    children: set[str] = set()

    if isinstance(node,list):
        for child in node:
            children.add(child)
    elif isinstance(node,dict):
        for child, descendants in node.items():
            children.add(child)
            if isinstance(descendants,list):
                children.update(descendants)
            elif isinstance(descendants,dict):
                children.update(get_all_children(child, node, visited))
    # for child in dicc.get(celltype, []):
    #     children |= get_all_children(child, dicc, visited.copy())
    
    return children

def save_anndata_files(adata_dicts, adata_dir):
    """
    Save per-sample AnnData objects as H5AD files.

    Parameters
    ----------
    adata_dicts : dict of str to anndata.AnnData
        Per-sample annotated data matrices.
    adata_dir : pathlib.Path
        Destination directory.

    Returns
    -------
    None
        Files are written as ``Sample_<sample>.h5ad``.
    """
    for k, adata in adata_dicts.items():
        adata.write(adata_dir / f"Sample_{k}.h5ad")

def labeling_cell_types(data_dicts, adata_dicts, cell_type_dict, sub_cell_type_dict, protein_markers, output_dir, custom_colors, save_anndata=True, arrows_orientation="contains"):
    """
    Assign primary and hierarchical phenotype labels to AnnData objects.

    Parameters
    ----------
    data_dicts : dict of str to pandas.DataFrame
        Cleaned cell-level tables used for marker comparisons.
    adata_dicts : dict of str to anndata.AnnData
        Corresponding per-sample AnnData objects.
    cell_type_dict : dict
        Primary marker-based phenotype definitions.
    sub_cell_type_dict : dict
        Hierarchical subtype definitions.
    protein_markers : list of str
        Primary protein-marker positivity columns.
    output_dir : pathlib.Path
        SpPrAn output directory.
    custom_colors : dict of str to str
        Phenotype color mapping.
    save_anndata : bool, default=True
        Whether updated AnnData objects should be written to disk.
    arrows_orientation : {"contains", "isContained"}, default="contains"
        Orientation used for hierarchy diagrams.

    Returns
    -------
    dict of str to anndata.AnnData
        AnnData objects containing primary labels, hierarchical subtype
        annotations, colors, and plotting metadata.

    Notes
    -----
    Primary population labels are stored in ``adata.obs["clusters"]``.
    Hierarchical subtype annotations are stored in additional observation
    columns associated with their parent phenotype.
    """
    list_final_subsets, list_cell_types, list_bigger_cell_types, custom_colors = sets_and_subsets_cell_types(cell_type_dict, sub_cell_type_dict, protein_markers, output_dir, custom_colors,arrows_orientation=arrows_orientation)

    for (k_data, data), (k, adata) in zip(data_dicts.items(), adata_dicts.items()):
        list_folder_names = ['General']
        new_obs_cols = {}
        # clusters = data.apply(assign_cell_type, axis=1, args=(cell_type_dict,)).tolist()
        clusters = assign_cell_types_vectorized(data, cell_type_dict, "General")
        new_obs_cols['clusters'] = clusters
        # dict_vl = {}
        for cell_type in cell_type_dict.keys():
            # new_obs_cols[f"only {cell_type}"] = [t if t==cell_type else "Other cells" for t in clusters]
            new_obs_cols[f"Only {cell_type}"] = clusters.where(clusters == cell_type, "Other cells")
        
            # if cell_type in list_final_subsets.keys():
            if cell_type in sub_cell_type_dict.keys():
                list_folder_names.append(cell_type)
                valid_labels = {cell_type} | get_all_children(cell_type, list_final_subsets)
                dict_to_pass = list_cell_types[cell_type]|{cell_type:list_cell_types['General'][cell_type]}
                protein_m_to_check = list(set([pm_to_check for ct_key,ct_rule in dict_to_pass.items() for pm_to_check in ct_rule.keys()]))
                # dict_vl[cell_type] = valid_labels
                # subtype_values = data.apply(assign_cell_type, axis=1, args=(list_cell_types[cell_type]|{cell_type:list_cell_types['General'][cell_type]},)).tolist()
                subtype_values = assign_cell_types_vectorized(data, dict_to_pass, cell_type)
                subtype_col = f"Subtypes of {cell_type}"
                new_obs_cols[subtype_col] = subtype_values
                for subcell_type in valid_labels:
                    if subcell_type==cell_type: continue
                    # new_obs_cols[f"only {subcell_type}"] = [t if t==subcell_type else "Other cells" for t in subtype_values]
                    # new_obs_cols[f"only {subcell_type}"] = subtype_values.where(subtype_values == subcell_type, "Other cells")
                    sct_temp = pd.Series("Other cells", index=subtype_values.index)
                    sct_temp = sct_temp.mask(subtype_values.isin(valid_labels),f"Other {cell_type}")
                    new_obs_cols[f"only {subcell_type}"] = sct_temp.mask(subtype_values==subcell_type,subtype_values)
                custom_colors[f"Other {cell_type}"] = custom_colors[cell_type]
        matches = {parent_key: list(set(child_dict.keys()) & set(sub_cell_type_dict.keys()))
            for parent_key, child_dict in sub_cell_type_dict.items()
            if isinstance(child_dict, dict)
        }
        matches = {k: v for k, v in matches.items() if v}
        ordered_keys = order_sub_cell_types(sub_cell_type_dict)
        for cell_type in ordered_keys:
            if cell_type not in cell_type_dict.keys():
                list_folder_names.append(cell_type)
                valid_labels = {cell_type} | get_all_children(cell_type, list_final_subsets)
                key_to_look = ""
                for key1,value1 in list_cell_types.items():
                    for key2 in value1.keys():
                        if key2==cell_type:
                            key_to_look = key1
                            break
                dict_to_pass = list_cell_types[cell_type]|{cell_type:list_cell_types[key_to_look][cell_type]}
                subtype_values = assign_cell_types_vectorized(data, dict_to_pass, cell_type)
                for type_of_cell in matches.keys():
                    if cell_type in matches[type_of_cell]:
                        if type_of_cell in cell_type_dict.keys():
                            temp_for_only = new_obs_cols[f"only {cell_type}"].mask(new_obs_cols[f"only {cell_type}"].ne(cell_type), "Other cells")
                        else:
                            temp_for_only = new_obs_cols[f"{type_of_cell} - only {cell_type}"].mask(new_obs_cols[f"{type_of_cell} - only {cell_type}"].ne(cell_type), "Other cells")
                new_obs_cols[f"Subtypes of {cell_type}"] = subtype_values.mask(temp_for_only.eq("Other cells"), "Other cells")
                new_obs_cols[f"Only {cell_type}"] = subtype_values.where(subtype_values == "Other cells", cell_type)
                for subcell_type in valid_labels:
                    if subcell_type==cell_type: continue
                    sct_temp = pd.Series("Other cells", index=subtype_values.index)
                    sct_temp = sct_temp.mask(subtype_values.isin(valid_labels),f"Other {cell_type}")
                    new_obs_cols[f"{cell_type} - only {subcell_type}"] = sct_temp.mask(subtype_values==subcell_type,subtype_values)
                custom_colors[f"Other {cell_type}"] = custom_colors[cell_type]
        # raise Exception (f"Debugging {k_data}:", print(adata.obs["clusters"].dtype), print(type(adata.obs["clusters"].iloc[0])), print(adata.obs["clusters"].head()), )
        new_obs = pd.DataFrame(new_obs_cols)
        new_obs.index = new_obs.index.map(str)
        adata.obs = new_obs if adata.obs.empty else adata.obs.join(new_obs)
        # for cell_type in final_subsets:
        #     valid_labels = {cell_type} | get_all_children(cell_type, final_subsets)
        #     adata.obs[f"{cell_type}"] = [t if t in valid_labels else "Other cells" for t in adata.obs['clusters']]
        # reverse_map = {t:bg for bg in bigger_cell_types for t in dict_vl[bg]}
        # adata.obs['bigger_cell_types'] = [t if t in bigger_cell_types else reverse_map.get(t) for t in adata.obs['clusters']]
        adata.uns['color_keys'] = [cell_type_name for cell_type_name in custom_colors.keys()]
        adata.uns['color_vals'] = [custom_colors[key] for key in adata.uns['color_keys']]
        adata.uns['folder_names'] = list_folder_names
        adata_dicts[k] = adata
    if save_anndata:
        adata_path = output_dir / "adata"
        adata_path.mkdir(parents=True, exist_ok=True)
        save_anndata_files(adata_dicts,adata_dir=adata_path)
    return adata_dicts

def create_or_load_anndata(config):
    """
    Load cached AnnData objects or generate them from source tables.

    Parameters
    ----------
    config : dict
        Validated SpPrAn configuration dictionary.

    Returns
    -------
    dict of str to anndata.AnnData
        Per-sample annotated data matrices.

    Notes
    -----
    Existing H5AD files are reused when allowed by the configuration. If files
    are missing, incomplete, or overwrite mode is enabled, SpPrAn rebuilds the
    AnnData objects from the original object-level input files.
    """
    output_dir = config['workspace']['output_dir']
    if not config['overwrite_existing_files']: adata_dicts = load_anndata_files(output_dir)
    else: adata_dicts = {}
    if any(v is None for v in adata_dicts.values()) or not adata_dicts:
        logger.warning("Generating anndata files for analysis...")
        file_names = filenames(config['workspace']['input_dir'], config['workspace']['filetype'])
        cell_type_dict = create_cell_type_dict(config['protein_markers'], config['cell_types'])
        sub_cell_type_dict = config['sub_cell_types'] if 'sub_cell_types' in config.keys() else {}
        data_dicts, adata_dicts = cleaned_data(
            file_names        = file_names, 
            output_dir        = output_dir, 
            protein_markers   = config['protein_markers'],
            cell_subtype_dict = sub_cell_type_dict,
            filetype          = config['workspace']['filetype'],
        )
        if len(adata_dicts)==0: raise ValueError("No anndata files were generated. Please review the configuration file.")
        logger.warning("Anndata generated!")
        adata_dicts = labeling_cell_types(data_dicts, adata_dicts, cell_type_dict, sub_cell_type_dict, config['protein_markers'], output_dir, config['custom_colors'], save_anndata = config['locally_save_anndata_files'], arrows_orientation=config['arrows_orientation'])
        logger.warning("Cells labeled by cell types/subtypes: DONE.")
        if config['locally_save_anndata_files']: logger.warning("Anndata saved!")
        else: logger.warning("Anndata not saved by user decision.")
    else: logger.warning("Anndata loaded!")
    return adata_dicts
