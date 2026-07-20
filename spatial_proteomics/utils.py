# Loading needed libraries

from pathlib import Path
import yaml
import warnings

import numpy as np
np.float_ = np.float64
import pandas as pd
from anndata import AnnData

import scanpy as sc
with warnings.catch_warnings(): # Filtering harmless warnings
    warnings.filterwarnings(
        "ignore",
        message=".*Dask DataFrame.*"
    )
    warnings.filterwarnings(
        "ignore",
        message=".*Import anndata.io.read_text instead.*"
    )
    import squidpy as sq

import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import ListedColormap, to_rgb, to_hex
import seaborn as sns

from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from scipy.stats import mannwhitneyu

from itertools import combinations, product
import networkx as nx

# Functions

## Loading configurations
def load_config(config_path):
    """
    Loads and validates a YAML configuration file.
    
    Parameters
    -----
    config_path : str | Path
        Path to the YAML configuration file.

    Return
    ------
    Dict
        Parsed configuration dictionary.
    """
    config_path = Path(config_path)
    if not config_path.exists(): raise FileNotFoundError(f"Config file not found: {config_path}")
    with open(config_path, "r") as file:
        try: config = yaml.safe_load(file)
        except yaml.YAMLError as e: raise ValueError(f"YAML parsing error in {config_path}") from e
    if not config: raise ValueError("Config file is empty")
    config['workspace']['input_dir'] = Path(config['workspace']['input_dir'])
    config['workspace']['output_dir'] = Path(config['workspace']['output_dir'])
    return config

## Saving adata files
def save_anndata_files(adata_dicts, adata_dir):
    """
    Saves AnnData files in anndata directory.
    
    Parameters
    -----
    adata_dicts : Dict[str:AnnData]
        Dictionary containing multiple annotated data matrices.
    adata_dir : Path
        Path to the anndata directory.
    """
    for k, adata in adata_dicts.items():
        adata.write(adata_dir / f"Sample_{k}.h5ad")

## Loading adata files
def load_anndata_files(output_dir):
    """
    Loads AnnData files from anndata directory (if they exists).
    
    Parameters
    -----
    output_dir : Path
        Path to the output directory.

    Return
    ------
    Dict[str:AnnData] | Dict[str:None] | Dict
       Dictionary containing multiple annotated data matrices. 
    """
    filepath = output_dir / "results" / "Samples Id.csv"
    if filepath.exists():
        sample_names = pd.read_csv(filepath)['Samples Id'].tolist()
        adata_dicts = {}
        for sample in sample_names:
            sample = str(sample)                                                                                             # Ensuring we can deal with numbers as identifiers.
            adata_path = output_dir / "adata" / f"Sample_{sample}.h5ad"
            if adata_path.exists():
                adata_dicts[sample] = sc.read_h5ad(adata_path)
            else:
                print(f"WARNING: Anndata loading cannot be completed. Sample_{sample}.h5ad doesn't exist or it's missing.") # An None value is saved in the dictionary.
                adata_dicts[sample] = None
        return adata_dicts
    print(f"WARNING: Anndata loading cannot be completed. File '{filepath}' doesn't exist or it's missing.")                # An empty dictionary is returned.
    return {}

## Obtaining filenames
def filenames(input_dir, filetype):
    """
    Gets files' paths from input directory.
    
    Parameters
    -----
    input_dir : Path
        Path to the input directory.
    filetype : str
        File extension to search.

    Return
    ------
    List[str]
       List containing multiple paths to input files. 
    """
    if input_dir.exists():
        entries = list(input_dir.iterdir())
        filenames = [entry for entry in entries if entry.is_file() and entry.match(f"*objects.{filetype}") and not entry.name.startswith("._") and not entry.name.startswith(".")]  # Exclude hidden files and macOS resource forks
        return filenames
    print(f"WARNING: Input directory '{input_dir}' doesn't exists or it's misspelled.") # An empty list is returned.
    return []

## Cleaning data
def cleaned_data(file_names, output_dir, filetype='tsv'):
    """
    Performs QC to input data and generates AnnData dictionary.
    
    Parameters
    -----
    file_names : List[str]
        List containing multiple paths to input files.
    output_dir : Path
        Path to the output directory.
    filetype : str
        File extension of all input files.

    Return
    ------
    data_dicts : Dict[str:pd.DataFrame]
       Dictionary containing all data cleaned. 
    adata_dicts : Dict[str:AnnData]
       Dictionary containing multiple annotated data matrices. 
    """
    adata_dicts = {}
    data_dicts  = {}
    if len(file_names)!=0:
        for filename in file_names:
            separator = '\t' if filetype=='tsv' else ','
            filename = str(filename)
            data = pd.read_csv(f"{filename}",sep=separator)
            if data.columns.isin(['Positivity - DAPI (MV - NUC)']).any():
                data = data[data['Positivity - DAPI (MV - NUC)']==1]                                                           # filter out cells without nucleus
            else: print("WARNING: This pipeline uses 'Positivity - DAPI' for filtering out cells without nucleus. However, such column was not found. Continue without filtering...")
            name_of_file = str(filename.split('/')[-1].split('\\')[-1])
            if name_of_file[:11]=='COMET_6x6 (':
                data['Name'] = [f"{s[11]}{s[14]}_{i:05}" for s,i in zip([name_of_file]*len(data),data.index)]   # retain letter and number for identifying samples
            else:
                data['Name'] = [f"{s}_{i:05}" for s,i in zip([name_of_file[:-12]]*len(data),data.index)]        # retain filename for identifying samples
            data.rename(columns={"Name":"cellID"}, inplace=True)
            
            pos_cols = data.columns[(data.columns.str.contains('Positivity'))&(data.columns.str.contains('MV'))]               # identify columns that says if protein marker is present (1) or abscent (0)
            data = data[~(data[pos_cols].eq(-1).any(axis=1))]                                                                  # filter out those cells that does contain a NaN value (-1) in the previous columns
            columns_nuc = data.columns[                                                                                        # keep columns with protein intensity in the nucleus
                (data.columns.isin(['cellID']))|
                ((data.columns.str.contains("MV - NUC - "))&(~data.columns.str.contains("Type")))
            ]
            columns_pos = data.columns[data.columns.isin(['cellID', 'X-coordinate', 'Y-coordinate'])]                          # keep spatial columns
            data_temp = data[columns_nuc].iloc[:,1:]
            adata = AnnData(                                                                                                   # generate AnnData file to include spatial data
                data_temp.set_index((str(x) for x in data_temp.index)),
                obsm={
                    "spatial": data[columns_pos].iloc[:,1:].to_numpy(),
                    "ID_cell":data[['cellID']].to_numpy(),
                    "Positivity":data[pos_cols].to_numpy()
                },
                uns={
                    "spatial":{"unique":{}}
                }
            )

            data = data[["cellID"]+pos_cols.tolist()]
            adata_dicts[f"{data.iloc[0,0][-8:-6]}"] = adata
            data_dicts[f"{data.iloc[0,0][-8:-6]}"] = data
        df = pd.DataFrame({"Samples Id": list(data_dicts.keys())})
        results_path = output_dir / "results"
        results_path.mkdir(parents=True, exist_ok=True)
        df.to_csv(results_path / "Samples Id.csv", index=False)
        
        df = pd.DataFrame({"Positivity column names": list(data[pos_cols].columns)})
        df.to_csv(results_path/ "Positivity column names.csv", index=False)
    else: 
        data_dicts = {}
        adata_dicts = {}
    return data_dicts, adata_dicts

## Considering groups and subgroups (its a reordering)
def ascending_order_for_cell_type_by_def(cell_type_dict):
    """
    Order the pre-established cell types considering types and subtypes of cells by definition.
    
    Parameters
    -----
    Dict[str:Dict]
       Dictionary containing the definition of each cell types.

    Return
    ------
    Dict[str:Dict]
       Dictionary containing the definition of each cell types ordered by cell subtypes.
    """
    
    return cell_type_dict

## Creating cell type dictionary
def create_cell_type_dict(protein_markers, cell_types):
    """
    Creates a dictionary with the definition of the pre-established cell types considering the presence or absent of protein markers.
    
    Parameters
    -----
    protein_markers : List[str]
        List containing protein markers.
    cell_types : Dict[str:List]
        Dictionary containing cell types as keys and a list with 0, 1, and None as values.

    Return
    ------
    Dict[str:Dict]
       Dictionary containing the definition of each cell types.
    """
    cell_type_dict = {}
    for cell_type, rule in cell_types.items():
        cell_type_dict[cell_type] = {protein_markers[i]:rule[i] for i in range(len(rule)) if rule[i] is not None}
    return cell_type_dict

def value_matches(observed, expected):
    """
    Check whether one observed marker value satisfies a rule.

    Rule meanings
    -------------
    0 or 1:
        The observed value must equal the expected value.
    list, tuple, or set:
        The observed value must be one of the allowed values.
    None:
        Ignore this marker.
    """
    if expected is None:
        return True
    if isinstance(expected, (list, tuple, set)):
        return observed in expected
    return observed == expected

## Assign cell types
def assign_cell_type(row, cell_type_dict):
    """
    Assigns a cell type to a cell.
    
    Parameters
    -----
    row : pd.Series()
        Row from DataFrame to modify.
    cell_type_dict : Dict[str:Dict]
        Dictionary containing the definition of each cell types.

    Return
    ------
    str
       Assigned cell type.
    """
    for cell_type, rule in cell_type_dict.items():
        if all(value_matches(row[marker], expected) for marker, expected in rule.items()):
            return str(cell_type)
    return "Other cells"

def assign_cell_types_vectorized(data, cell_type_dict, container):
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

def get_index_nones(cell_types, protein_markers, sub_protein_markers, sub_cell_types):
    """
    Identifies the indices of None values in the cell types definition.

    Parameters
    -----
    cell_types : Dict[str:List]
        Dictionary containing cell types as keys and a list with 0, 1, and None as values.
    protein_markers : List[str]
        List containing protein markers.
    sub_protein_markers : List[str]
        List containing sub-protein markers.
    sub_cell_types : Dict[str:Dict]
        Dictionary containing the definition of each sub cell types.

    Return
    ------
    Dict[str:List]
       Dictionary containing the indices of None values for each cell type.
    """
    none_dicc = {}
    total_markers = len(protein_markers)
    none_subindices = [index+total_markers for index in range(len(sub_protein_markers['all']))]
    for cell_type in cell_types.keys():
        none_indices = [index for index, value in enumerate(protein_markers) if value not in cell_types[cell_type].keys()]
        none_dicc[cell_type]=none_indices
    for cell_type in sub_cell_types.keys():
        none_indices_temp = none_dicc[cell_type]+none_subindices
        for subcell_type in sub_cell_types[cell_type].keys():
            # none_indices = [index for index in none_indices_temp if index not in [total_markers + sub_protein_markers['all'].index(val) for val in sub_cell_types[cell_type][subcell_type].keys()]]
            none_indices = [total_markers + sub_protein_markers['all'].index(val) for val in sub_protein_markers[cell_type] if val not in sub_cell_types[cell_type][subcell_type].keys()]
            none_dicc[subcell_type]=none_indices
    return none_dicc

def fill_n_order_df_by_total_none_values(df, idx_name, order_dicc):
    """
    Orders a DataFrame based on the number of None values in each column.

    Parameters
    -----
    df : pd.DataFrame
        DataFrame to order.
    idx_name : str
        Name of the index to order by.
    order_dicc : Dict[str:List]
        Dictionary containing the indices of None values for each cell type.

    Return
    ------
    pd.DataFrame
       Ordered DataFrame.
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

def expand_subsets(col,df):
    """
    Expands the subsets of a cell type based on the presence or absence of protein markers.

    Parameters
    -----
    col : str
        Column name in the DataFrame.
    df : pd.DataFrame
        DataFrame containing the cell type information.

    Return
    ------
    pd.DataFrame
       Expanded DataFrame with subsets of the cell type.
    """
    n = df[col].apply(lambda x: isinstance(x, list)).sum()
    combos = list(product(*([ [0,1] ] * n)))
    list_to_df = []
    aux = 0
    for val in df[col]:
        if isinstance(val,list):
            list_to_df.append([c[aux] for c in combos])
            aux+=1
        else: list_to_df.append([val]*2**n)
    return pd.DataFrame(list_to_df, columns=[f"{col}_{''.join(map(str, c))}" for c in combos])

def move_bigger_groups(dic, zero_degree_items, target, final_subsets):
    """
    Moves bigger groups to the target level in the dictionary.

    Parameters
    -----
    dic : Dict[int:List]
        Dictionary containing groups of cell types.
    zero_degree_items : List[str]
        List of cell types with zero degrees.
    target : int
        Target level in the dictionary.
    final_subsets : Dict[str:List] 
        Dictionary containing the final subsets of cell types.  
    
    Return
    ------
    Dict[int:List]
       Dictionary containing the updated groups of cell types.
    int
         Maximum level in the updated dictionary.
    """

    find = []
    for k, item_list in dic.items():
        for item in zero_degree_items:
            if item in item_list:
                item_list.remove(item)
                find.append(item)
    dic[target].extend(find)
    
    # Number of None values will determine the order of the cell types. The more None values, the bigger the group.
    dicc_target = {k:[] for k in reversed(range(target+1))}
    dicc_target[target] = dic[target]
    for i in reversed(range(target)):
        # if dic[1]!=0:
        for item in dic[i]:
            check = 0
            for k, item_list in final_subsets.items():
                if item in item_list:
                    for keys, vals in dicc_target.items():
                        if k in vals:
                            check+=1
                            dicc_target[keys-1].append(item)
                            break
            if check == 0: dicc_target[target].append(item)
    for keys, vals in dicc_target.items():
        dicc_target[keys] = list(dict.fromkeys(vals)) # Remove duplicates while preserving order
        for val in vals:
            if val in final_subsets.keys():
                dicc_target[keys-1].extend(final_subsets[val])
    dicc_target = {k:v for k,v in dicc_target.items() if len(v)!=0}
    if len(dicc_target)<target+1:
        for i in range(target+1):
            if i not in dicc_target:
                for j in range(i+1, target+1):
                    if j in dicc_target:
                        dicc_target[i] = dicc_target[j]
                        dicc_target.pop(j)
                        break
    return dicc_target, max([k for k in dicc_target.keys()])

def mix_hex_colors(colors,dim_factor=0.25,background="#ffffff"):
    """
    Combine multiple HEX colors and optionally dim/desaturate the result
    by blending it toward a background color.

    Parameters
    ----------
    colors : list[str]
        HEX colors to combine.
    dim_factor : float
        Amount of blending toward the background.
        0.0 keeps the original mixed color.
        1.0 becomes the background color.
    background : str
        Background color used for dimming.
        Use '#ffffff' for white plots or '#000000' for dark plots.

    Returns
    -------
    str
        Combined HEX color.
    """
    rgb_colors = [to_rgb(color) for color in colors]
    mixed_rgb = tuple(sum(color[channel] for color in rgb_colors) / len(rgb_colors) for channel in range(3))
    background_rgb = to_rgb(background)
    dimmed_rgb = tuple((1 - dim_factor) * mixed + dim_factor * background_value for mixed, background_value in zip(mixed_rgb, background_rgb))
    return to_hex(dimmed_rgb)

def create_intersection_color(intersection_label, list_intersected_cell_subtypes, custom_colors, dim_factor=0.20, background="#ffffff"):
    """
    Generate a color for an intersection label from its component colors.

    Parameters
    ----------
    intersection_label : str
        Cell type that comes from the intersection of two cell types.
    list_intersected_cell_subtypes: List[str]
        List with the cell types who are intersected.
    custom_colors : Dict[str:str]
       Dictionary containing HEX colors per cell type. 
    dim_factor : float
        Amount of blending toward the background.
        0.0 keeps the original mixed color.
        1.0 becomes the background color.
    background : str
        Background color used for dimming.
        Use '#ffffff' for white plots or '#000000' for dark plots.

    Returns
    -------
    str
        Combined HEX color.
    """
    if intersection_label not in custom_colors.keys():
        component_colors = [custom_colors[subtype] for subtype in list_intersected_cell_subtypes]
        return mix_hex_colors(component_colors, dim_factor=dim_factor, background=background)
    else: return custom_colors[intersection_label]

def plot_digraph(ct_df, protein_markers, output_dir, output_filename_identifier, order_dicc, custom_colors):
    """
    Plots a directed graph of the final subsets of cell types.

    Parameters
    -----
    ct_df : pd.DataFrame
        DataFrame containing the cell type data.
    protein_markers : List[str]
        List of protein markers.
    output_dir : Path
        Path to the output directory.
    output_filename_identifier : str
        Identifier for the output filename.
    order_dicc : Dict[int, List]
        Ordered dictionary of cell types.
    custom_colors : Dict[str, str]
        Dictionary mapping cell types to colors.

    Return
    -----
    custom_colors : Dict[str, str]
        Updated dictionary mapping cell types to colors.
    """
    combinations_list = [list(c) for c in combinations(ct_df.columns, 2)]
    dicc_subsets = {col:expand_subsets(col, ct_df) for col in ct_df.columns}
    
    dicc_temp = {}
    max_N = max([k for k in order_dicc.keys()])
    for total_N, cell_types in order_dicc.items():
        dicc_temp[total_N] = []
        for cell_type in cell_types:
            if cell_type in ct_df.columns:
                dicc_temp[total_N].append(cell_type)
        if total_N==max_N and len(dicc_temp[total_N])==0:
            dicc_temp.pop(total_N)
            max_N = max([k for k in order_dicc.keys() if k<total_N])

    cell_type_chart = nx.DiGraph()
    for cell_type in ct_df.columns: cell_type_chart.add_node(cell_type) 
    
    final_subsets = {cell_type:[] for cell_type in ct_df.columns}
    for elem in combinations_list:
        temp_0 = list(dicc_subsets[elem[0]].columns)
        temp_1 = [col for col in dicc_subsets[elem[0]].columns for col_c in dicc_subsets[elem[1]].columns if dicc_subsets[elem[0]][col].equals(dicc_subsets[elem[1]][col_c])]
        match len(temp_1):
            case x if x==len(temp_0): 
                print(f'{elem[0]} is subset of {elem[1]}')
                cell_type_chart.add_edges_from([(elem[0],elem[1])])
                # print(f"We added the path from {elem[0]} to {elem[1]}") # debugging
                final_subsets[elem[1]].append(elem[0])
                if len(temp_1) != 1:
                    for pos_subset in range([key for key in dicc_temp.keys() if elem[0] in dicc_temp[key]][0]):
                        for pot_subset in dicc_temp[pos_subset]:
                            if nx.has_path(cell_type_chart,source=pot_subset,target=elem[1]) and nx.shortest_path_length(cell_type_chart, source=pot_subset, target=elem[1])==1:
                                # print(f"We attempted removing the path from {pot_subset} to {elem[1]}") # debugging
                                cell_type_chart.remove_edge(pot_subset, elem[1])
                                # print(f"We succeeded removing the path from {pot_subset} to {elem[1]}") # debugging
                                final_subsets[elem[1]].remove(pot_subset)
            case x if x>0: 
                whos_who = []
                for pos_subset in range([key for key in dicc_temp.keys() if elem[0] in dicc_temp[key]][0]):
                    for pot_subset in dicc_temp[pos_subset]:
                        if nx.has_path(cell_type_chart,source=pot_subset,target=elem[0]) and nx.has_path(cell_type_chart,source=pot_subset,target=elem[1]):
                            whos_who.append(pot_subset)
                list_temp_temp = list(set(col for col in temp_1 for col_s in whos_who for col_c in dicc_subsets[col_s].columns if dicc_subsets[elem[0]][col].equals(dicc_subsets[col_s][col_c])))
                temp_new = [x for x in temp_1 if x not in list_temp_temp]
                limit_dim_color = {k:0.2+0.05*k for k in range(len(temp_new))} if len(temp_new) <= 14 else {k:0.2+(0.7*k/len(temp_new)) for k in range(len(temp_new))}
                0.2+len(temp_new)
                match len(temp_new):
                    case x if x>1:
                        df_temp = dicc_subsets[elem[0]].rename(index={idx:val for idx,val in enumerate(protein_markers)})
                        print(f"There are sub cell types that are not defined by user and belong to both {elem[0]} and {elem[1]}.")#\nThey are the following:")
                        
                        for sct_def in temp_new:
                            # print(df_temp[sct_def].to_dict())
                            if sct_def in cell_type_chart:
                                if not nx.has_path(cell_type_chart,source=sct_def,target=elem[0]): 
                                    cell_type_chart.add_edges_from([(sct_def,elem[0])])
                                    final_subsets[elem[0]].append(sct_def)
                                if not nx.has_path(cell_type_chart,source=sct_def,target=elem[1]): 
                                    cell_type_chart.add_edges_from([(sct_def,elem[1])])
                                    final_subsets[elem[1]].append(sct_def)
                            else: 
                                cell_type_chart.add_node(sct_def) # Adding the node to the graph if it doesn't exist # debugging
                                cell_type_chart.add_edges_from([(sct_def,elem[0])])
                                cell_type_chart.add_edges_from([(sct_def,elem[1])])
                                final_subsets[elem[0]].append(sct_def)
                                final_subsets[elem[1]].append(sct_def)
                            custom_colors[sct_def] = create_intersection_color(sct_def,[elem[0],elem[1]],custom_colors,limit_dim_color[temp_new.index(sct_def)])
                    case x if x==1:
                        df_temp = dicc_subsets[elem[0]].rename(index={idx:val for idx,val in enumerate(protein_markers)})
                        print(f"There is a sub cell type that is not defined by user and belongs to both {elem[0]} and {elem[1]}.")#\nIt is the following:")
                        # print(df_temp[temp_new[0]].to_dict())
                        if temp_new[0] in cell_type_chart:
                            if not nx.has_path(cell_type_chart,source=temp_new[0],target=elem[0]): 
                                cell_type_chart.add_edges_from([(temp_new[0],elem[0])])
                                final_subsets[elem[0]].append(temp_new[0])
                            if not nx.has_path(cell_type_chart,source=temp_new[0],target=elem[1]): 
                                cell_type_chart.add_edges_from([(temp_new[0],elem[1])])
                                final_subsets[elem[1]].append(temp_new[0])
                        else:
                            cell_type_chart.add_node(temp_new[0]) # Adding the node to the graph if it doesn't exist # debugging  
                            cell_type_chart.add_edges_from([(temp_new[0],elem[0])])
                            cell_type_chart.add_edges_from([(temp_new[0],elem[1])])
                            final_subsets[elem[0]].append(temp_new[0])
                            final_subsets[elem[1]].append(temp_new[0])
                        custom_colors[temp_new[0]] = create_intersection_color(temp_new[0],[elem[0],elem[1]],custom_colors)
                    case _: pass
            case _: pass
    
    possible_duplicates = {cell_type:list_val for cell_type, list_val in final_subsets.items() if len(list_val)>1}
    for cell_type, list_val in possible_duplicates.items():
        for ct in list_val:
            for ct_temp in [k for k in possible_duplicates.keys() if k!=cell_type]:
                if set([cell_type,ct]).issubset(set(possible_duplicates[ct_temp])):
                    final_subsets[ct_temp].remove(ct)
                    cell_type_chart.remove_edge(ct, ct_temp)
    csv_path = output_dir / "results" / output_filename_identifier
    if not csv_path.exists(): csv_path.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({k:[v] for k,v in final_subsets.items()},index=['Subsets']).T.to_csv(csv_path / "Cell type subsets.csv", index=True)
    
    # ct_final_subsets = pd.DataFrame(cell_types, index=protein_markers, dtype=object)
    ct_final_subsets = ct_df.copy()
    final_order = []
    new_columns = {}
    cols_ct_final_subsets = list(ct_final_subsets.columns)
    for k,val in final_subsets.items():
        keys_new_columns = list(new_columns.keys())
        columns_to_check = cols_ct_final_subsets+keys_new_columns
        if val:
            for ct in val:
                if ct not in columns_to_check:
                    new_columns[ct]=dicc_subsets[k][ct].rename(index={idx:pmarker for idx,pmarker in enumerate(protein_markers)})
                    final_order.append(ct)
        final_order.append(k)
    if new_columns: ct_final_subsets = pd.concat([ct_final_subsets,pd.DataFrame(new_columns)],axis=1)
    ct_final_subsets = ct_final_subsets[final_order]
    # ct_final_subsets = {col:ct_final_subsets[col].to_list() for col in ct_final_subsets.columns}
    ct_final_subsets = ct_final_subsets.to_dict()
    # print("ct_final_subsets\n",ct_final_subsets)
    df_df_temp = pd.DataFrame(ct_final_subsets, dtype=object)
    none_mask_temp = df_df_temp.isna()
    for col in df_df_temp.columns:
        if none_mask_temp[col].any():
            df_df_temp.loc[df_df_temp[col].isna(), col] = [[0,1] for _ in range(df_df_temp[col].isna().sum())]
    df_df_temp.to_csv(csv_path / "All cell type definitions.csv", index=True)
    
    # plt.figure(figsize=(6,6),dpi=300)
    # cell_type_chart = cell_type_chart.subgraph(list(ct_df.columns)).copy() # This line avoids to plot not-defined cell types
    # pos = nx.spring_layout(cell_type_chart, seed=538610)  # This can be replaced with other layout strategies
    pos = {}
    zero_degree_ct = [cell_type for cell_type in ct_df.columns if cell_type_chart.degree[cell_type] == 0]
    if max_N != 0: dicc_temp, max_N = move_bigger_groups(dicc_temp, zero_degree_ct, max_N, final_subsets)
    max_per_group = max([len(N_pos) for N_pos in dicc_temp.values()])
    for i in reversed(range(max_N+1)):
        for j in reversed(range(len(dicc_temp[i]))):
            cell_type = dicc_temp[i][j]
            # pos[cell_type] = np.array([j+(max_per_group/len(dicc_temp[i])),i+(max_per_group/len(dicc_temp))]) if len(dicc_temp[i])!=1 else np.array([j+((max_per_group+1)/2),i+(max_per_group/len(dicc_temp))])
            pos[cell_type] = (j+((max_per_group-len(dicc_temp[i]))/2),i)
    nx.draw(cell_type_chart, pos, with_labels=True,
            node_size=8000,
            labels = {node:node.replace(" ", "\n") for node in cell_type_chart.nodes()},
            node_color=[custom_colors[node] if node in custom_colors.keys() else 'lightgray' for node in cell_type_chart.nodes()],
            font_size=16,
            width = 3,
            margins = 0.15,
            edge_color='gray', arrowstyle='-|>')
    plt.title("Cell types hierarchy", size=20)
    plots_path = output_dir / "plots" / output_filename_identifier 
    if not plots_path.exists(): plots_path.mkdir(parents=True, exist_ok=True)
    plt.savefig(plots_path / "Cell types hierarchy.png")
    plt.close()
    return final_subsets, ct_final_subsets, dicc_temp[max_N], custom_colors
    
def sets_and_subsets_cell_types(cell_type_dict, sub_cell_type_dict, protein_markers, output_dir, custom_colors):
    """
    Identifies set and subsets of cell types based on their marker expressions.
    
    Parameters
    -----
    cell_type_dict : Dict[str:Dict]
        Dictionary containing the definition of each cell types.
    sub_cell_type_dict : Dict[str:Dict]
        Dictionary containing the definition of each sub cell types.
    protein_markers : List[str]
        List containing protein markers.

    Return
    ------
    final_subsets : Dict[str:List]
        Dictionary containing the final subsets of cell types.
    ct_final_subsets : Dict[str:Dict]
        Dictionary containing the final definition of each cell types.
    dicc_temp[max_N] : Lists[str]
        List with the biggest cell types in the analysis.
    custom_colors : Dict[str, str]
        Updated dictionary mapping cell types to colors.
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
        full_index = protein_markers + [p for p in sub_protein_markers[sct] if p not in protein_markers]
        df = pd.DataFrame(data, index=full_index, dtype=object)
        for prot in protein_markers: df.loc[prot] = cell_type_dict[sct].get(prot, None)
        df, df_og = fill_n_order_df_by_total_none_values(df, 'Protein Marker', dicc_temp)
        sct_dfs[sct] = df
        # sct_dfs_og[sct] = df_og

    list_final_subsets = {}
    list_sub_cell_types = {}
    list_bigger_cell_types = {}
    list_final_subsets["General"], list_sub_cell_types["General"], list_bigger_cell_types["General"], custom_colors = plot_digraph(ct_df, protein_markers, output_dir, "General", dicc_temp, custom_colors) 
    for sct, df in sct_dfs.items():
        list_final_subsets[sct], list_sub_cell_types[sct], list_bigger_cell_types[sct], custom_colors = plot_digraph(df, list(df.index), output_dir, sct, dicc_temp, custom_colors)
    return list_final_subsets, list_sub_cell_types, list_bigger_cell_types, custom_colors
    # return final_subsets, ct_final_subsets, dicc_temp[max_N]

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

    Return
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

## Labeling cell types
def labeling_cell_types(data_dicts, adata_dicts, cell_type_dict, sub_cell_type_dict, protein_markers, output_dir, custom_colors, save_anndata=True):
    """
    Labels all cells in AnnData with the corresponding cell type.
    
    Parameters
    -----
    data_dicts : Dict[str:pd.DataFrame]
       Dictionary containing all data cleaned. 
    adata_dicts : Dict[str:AnnData]
       Dictionary containing multiple annotated data matrices. 
    cell_type_dict : Dict[str:Dict]
        Dictionary containing the definition of each cell types.
    sub_cell_type_dict : Dict[str:Dict]
        Dictionary containing the definition of each sub cell types.
    output_dir : Path
        Path to the output directory.
    save_anndata : Bool (Optional; Default is True)
        Boolean value to decide if AnnData files will be saved in anndata directory or not.

    Return
    ------
    adata_dicts : Dict[str:AnnData]
       Dictionary containing multiple updated annotated data matrices.
    """
    list_final_subsets, list_cell_types, list_bigger_cell_types, custom_colors = sets_and_subsets_cell_types(cell_type_dict, sub_cell_type_dict, protein_markers, output_dir, custom_colors)

    for (k_data, data), (k, adata) in zip(data_dicts.items(), adata_dicts.items()):
        list_folder_names = ['General']
        new_obs_cols = {}
        # clusters = data.apply(assign_cell_type, axis=1, args=(cell_type_dict,)).tolist()
        clusters = assign_cell_types_vectorized(data, cell_type_dict, "General")
        new_obs_cols['clusters'] = clusters
        # dict_vl = {}
        for cell_type in cell_type_dict.keys():
            # new_obs_cols[f"only {cell_type}"] = [t if t==cell_type else "Other cells" for t in clusters]
            new_obs_cols[f"only {cell_type}"] = clusters.where(clusters == cell_type, "Other cells")
        
            if cell_type in list_final_subsets.keys():
                if cell_type in sub_cell_type_dict.keys():
                    list_folder_names.append(cell_type)
                    valid_labels = {cell_type} | get_all_children(cell_type, list_final_subsets)
                    # dict_vl[cell_type] = valid_labels
                    # subtype_values = data.apply(assign_cell_type, axis=1, args=(list_cell_types[cell_type]|{cell_type:list_cell_types['General'][cell_type]},)).tolist()
                    subtype_values = assign_cell_types_vectorized(data, list_cell_types[cell_type]|{cell_type:list_cell_types['General'][cell_type]},cell_type)
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

## Create or load anndata
def create_or_load_anndata(config):
    """
    Creates or loads AnnData files from anndata directory (if they exists).
    
    Parameters
    -----
    config : Dict
        Parsed configuration dictionary.

    Return
    ------
    adata_dicts : Dict[str:AnnData]
       Dictionary containing multiple annotated data matrices. 
    """
    output_dir = config['workspace']['output_dir']
    if not config['overwrite_existing_files']: adata_dicts = load_anndata_files(output_dir)
    else: adata_dicts = {}
    if any(v is None for v in adata_dicts.values()) or not adata_dicts:
        print("Generating anndata files for analysis...")
        file_names = filenames(config['workspace']['input_dir'], config['workspace']['filetype'])
        data_dicts, adata_dicts = cleaned_data(
            file_names   = file_names, 
            output_dir   = output_dir, 
            filetype     = config['workspace']['filetype']
        )
        if len(adata_dicts)==0: raise ValueError("No anndata files were generated. Please check your input directory and configuration file.")
        print("Anndata generated!")
        cell_type_dict = create_cell_type_dict(config['protein_markers'], config['cell_types'])
        sub_cell_type_dict = config['sub_cell_types'] if 'sub_cell_types' in config.keys() else {}
        adata_dicts = labeling_cell_types(data_dicts, adata_dicts, cell_type_dict, sub_cell_type_dict, config['protein_markers'], output_dir, config['custom_colors'], save_anndata = config['locally_save_anndata_files'])
        if config['locally_save_anndata_files']: print("Anndata saved!")
        else: print("WARNING: Anndata not saved.")
    else: print("Anndata loaded!")
    return adata_dicts


def plotting_helper(custom_colors, overwrite_existing_files, dpi, size, plots_path, k, adata, folder_name):
    if folder_name == 'General':
        universe_to_color = 'clusters' 
        title_name = f"Sample {k} ({adata.n_obs} cells)"
    else:
        universe_to_color = f"Subtypes of {folder_name}"
        total_subcells = (adata.obs['clusters'] == folder_name).sum()
        title_name = f"Sample {k} - {folder_name} ({total_subcells} cells)"
    save_namefile = plots_path / f"Spatial - {title_name}.png"

    if save_namefile.exists() and not overwrite_existing_files:
        print(f"File 'Spatial - {title_name}.png' already exists...")
        for cell_type, selected_color in custom_colors.items():
            if cell_type=='Other cells': continue
            if cell_type.startswith('Other'): continue
            if cell_type not in adata.obs[universe_to_color].values: continue
            if (adata.obs[universe_to_color] == cell_type).sum()==0: continue
            title_name2 = f"Sample {k} - {cell_type} ({adata.n_obs} cells)" if folder_name == 'General' else f"Sample {k} - {folder_name} - {cell_type} ({total_subcells} cells)"
            save_namefile2 = plots_path / f"Spatial - {title_name2}.png"
            if save_namefile2.exists() and not overwrite_existing_files:
                print(f"File 'Spatial - {title_name2}.png' already exists...")
                continue
            # adata.obs[f"only {cell_type}"]
            # adata.obs[f"Subtypes of {cell_type}"] # Sub universe
            color_spatial = f"only {cell_type}"
                # cmap_gene = selected_color
            # own_palette_list = [selected_color,custom_colors['Other cells']] if cell_type<'Other cells' else [custom_colors['Other cells'],selected_color]
            own_palette_list = [custom_colors[label] for label in sorted([cell_type, f"Other cells"])] if folder_name=='General' else [custom_colors[label] for label in sorted([cell_type, f"Other {folder_name}", f"Other cells"])]
            own_palette = ListedColormap(own_palette_list)
                # cmap_gene = ListedColormap(cmap_gene)
                
            fig, ax = plt.subplots()
            with warnings.catch_warnings(): # Filtering harmless warning
                warnings.filterwarnings(
                        "ignore",
                        message=".*No data for colormapping provided via 'c'. Parameters 'cmap', 'norm' will be ignored.*"
                    )
                sq.pl.spatial_scatter(
                        adata, 
                        shape=None, 
                        library_id="unique",
                        color=color_spatial,
                        title=title_name2,
                        dpi=dpi,
                        # cmap=cmap_gene,
                        palette=own_palette,
                        size=size,
                        ax=ax)
            ax.set_facecolor("black")
            ax.set_aspect('auto')
            # ax.set_aspect('equal')
            fig.subplots_adjust(right=0.8)
            # fig.tight_layout()
            plt.savefig(save_namefile2, dpi=dpi)
            plt.close(fig)
            # print(f"File 'Spatial - {title_name2}.png' created!")
        # print(f"Spatial plots for the cell subtypes of each single cell types in {folder_name} has been created!")
        return
    color_spatial = universe_to_color
        # cmap_gene = None
    which_colors = []
    colors_temp = list(custom_colors.keys())
    colors_temp.sort()
    custom_colors_t = {k:custom_colors[k] for k in colors_temp}
    for cell_type,selected_color in custom_colors_t.items():
        if (adata.obs[universe_to_color] == cell_type).sum()!=0: which_colors.append(selected_color)
    own_palette_list = (which_colors)
    own_palette = ListedColormap(own_palette_list)
        
    fig, ax = plt.subplots()
    with warnings.catch_warnings(): # Filtering harmless warning
        warnings.filterwarnings(
                "ignore",
                message=".*No data for colormapping provided via 'c'. Parameters 'cmap', 'norm' will be ignored.*"
            )
        sq.pl.spatial_scatter(
                adata, 
                shape=None, 
                library_id="unique",
                color=color_spatial,
                title=title_name,
                dpi=dpi,
                # cmap=cmap_gene,
                palette=own_palette,
                size=size,
                ax=ax)
    ax.set_facecolor("black")
    ax.set_aspect('auto')
    # ax.set_aspect('equal')
    fig.subplots_adjust(right=0.8)
    # plt.legend(ncol=1,loc='center right',bbox_to_anchor=(0.5, -0.05))
    # fig.tight_layout()
    plt.savefig(save_namefile, dpi=dpi)
    plt.close(fig)
    # print(f"File 'Spatial - {title_name}.png' created!")
    # print(f"Spatial plots for the cell types in {folder_name} has been created!")
        
    for cell_type, selected_color in custom_colors.items():
        if cell_type=='Other cells': continue
        if cell_type.startswith('Other'): continue
        if cell_type not in adata.obs[universe_to_color].values: continue
        if (adata.obs[universe_to_color] == cell_type).sum()==0: continue
            
        title_name2 = f"Sample {k} - {cell_type} ({adata.n_obs} cells)" if folder_name == 'General' else f"Sample {k} - {folder_name} - {cell_type} ({total_subcells} cells)"
        save_namefile2 = plots_path / f"Spatial - {title_name2}.png"
        if save_namefile2.exists() and not overwrite_existing_files: 
            print(f"File 'Spatial - {title_name2}.png' already exists...")
            continue
        color_spatial = f"only {cell_type}"
            # cmap_gene = selected_color
        # own_palette_list = [selected_color,custom_colors['Other cells']] if cell_type<'Other cells' else [custom_colors['Other cells'],selected_color]
        own_palette_list = [custom_colors[label] for label in sorted([cell_type, f"Other cells"])] if folder_name=='General' else [custom_colors[label] for label in sorted([cell_type, f"Other {folder_name}", f"Other cells"])]
        own_palette = ListedColormap(own_palette_list)
            # cmap_gene = ListedColormap(cmap_gene)
        fig, ax = plt.subplots()
        with warnings.catch_warnings(): # Filtering harmless warning
            warnings.filterwarnings(
                    "ignore",
                    message=".*No data for colormapping provided via 'c'. Parameters 'cmap', 'norm' will be ignored.*"
                )
            sq.pl.spatial_scatter(
                    adata, 
                    shape=None, 
                    library_id="unique",
                    color=color_spatial,
                    title=title_name2,
                    dpi=dpi,
                    # cmap=cmap_gene,
                    palette=own_palette,
                    size=size,
                    ax=ax)
        ax.set_facecolor("black")
        ax.set_aspect('auto')
        # ax.set_aspect('equal')
        fig.subplots_adjust(right=0.8)
        # fig.tight_layout()
        plt.savefig(save_namefile2, dpi=dpi)
        plt.close(fig)
        # print(f"File 'Spatial - {title_name2}.png' created!")
    # print(f"Spatial plots for the cell subtypes of each single cell types in {folder_name} has been created!")

## Plotting spacial data
def plot_spatial(adata_dicts,custom_colors,output_dir,overwrite_existing_files=False,dpi=300,size=50):
    """
    Plots spatial data from AnnData and saves it in plots directory.
    
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
    dpi : int (Optional; Default is 300)
        Dots per inch for spatial plot.
    size : int (Optional; Default is 50)
        Scatter dots size for spatial plot.
    """
    rcParams["figure.figsize"] = (10,10)
    plots_path = output_dir / "plots" 
    plots_path.mkdir(parents=True, exist_ok=True)
    for k, adata in adata_dicts.items():
        for key,val in zip(adata.uns['color_keys'], adata.uns['color_vals']):
            if key not in custom_colors.keys(): custom_colors[key] = val
        for folder_name in adata.uns['folder_names']:
            plots_path_folder = plots_path / folder_name
            plots_path_folder.mkdir(parents=True, exist_ok=True)
            plotting_helper(custom_colors, overwrite_existing_files, dpi, size, plots_path_folder, k, adata, folder_name)
        print(f"All spatial plots for cell subtypes in Sample {k} have been created!")
        # raise Exception (f"NotFinishedFunctionError: Just testing. This is custom_colors: {custom_colors}")

## Calculating cell proportions and saving them in csv files
def calculate_cell_proportions(adata_dicts,custom_colors, output_dir, overwrite_existing_files=False):
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
    for k, adata in adata_dicts.items():
        for key,val in zip(adata.uns['color_keys'], adata.uns['color_vals']):
            if key not in custom_colors.keys(): custom_colors[key] = val
        for folder_name in adata.uns['folder_names']:
            if folder_name == 'General':
                universe_to_count = adata.obs['clusters'] 
                total_subcells = adata.n_obs
                title_name = f"Sample {k} ({adata.n_obs} cells)"
            else:
                universe_to_count = adata.obs[f"Subtypes of {folder_name}"][adata.obs[f"Subtypes of {folder_name}"]!="Other cells"]
                total_subcells = (adata.obs['clusters'] == folder_name).sum()
                title_name = f"Sample {k} - {folder_name} ({total_subcells} cells)"
            save_namefile = output_dir / "results" /  folder_name / f"Cell type proportions - {title_name}.csv"
            if save_namefile.exists() and not overwrite_existing_files: 
                print(f"File 'Cell type proportions - {title_name}' already exists...")
                continue
            sum_dict = {}
            sum_dict["Total cells"] = adata.n_obs
            if sum_dict["Total cells"]==0: 
                print(f"No cells in '{title_name}' file. Skipping...")
                break
            prop_dict = {}
            for cell_type,_ in custom_colors.items():
                if cell_type not in universe_to_count.values: continue
                sum_dict[cell_type]  = (universe_to_count == cell_type).sum()
                prop_dict[cell_type] = str(np.round(sum_dict[cell_type]/total_subcells*100,2))+"%" if total_subcells != 0 else "0%"
            
            df = pd.DataFrame(data = {cell_type:[sum_dict[cell_type],prop_dict[cell_type]] for cell_type,_ in custom_colors.items() if cell_type in universe_to_count.values}, index = ["Total cells in cell type", "Percentage"])
            df.to_csv(save_namefile, index=False)
        print(f"All 'Cell type proportions' csv files for Sample {k} have been created!")









