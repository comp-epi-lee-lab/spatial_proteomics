# Loading needed libraries

from pathlib import Path
import yaml

import numpy as np
np.float_ = np.float64
import pandas as pd
from anndata import AnnData

import scanpy as sc
import squidpy as sq

import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import ListedColormap
import seaborn as sns

from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from scipy.stats import mannwhitneyu

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
            adata_path = output_dir / "adata" / f"Sample_{sample}.h5ad"
            if adata_path.exists():
                adata_dicts[sample] = sc.read_h5ad(adata_path)
            else:
                print(f"Unexpected error: Sample_{sample}.h5ad doesn't exists or it's missing.")
                adata_dicts[sample] = None
        return adata_dicts
    print(f"File '{filepath}' doesn't exists or it's missing. An empty dictionary is returned.")
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
        filenames = [entry for entry in entries if entry.is_file() and entry.match(f"*objects.{filetype}")]
        return filenames
    print("Input directory doesn't exists or it's missing. An empty list is returned.")
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
            data = pd.read_csv(f"{filename}",sep=separator)
            data = data[data['Positivity - DAPI (MV - NUC)']==1]                                        # filter out cells without nucleus
            
            temp = pd.concat([data.iloc[:,3], data.iloc[:, 12:]], axis=1)                               # retain only "Name" and data columns
            temp['Name'] = [f"{s[11]}{s[14]}_{i:05}" for s,i in zip(temp['Name'],temp['Name'].index)]   # retain letter and number for identifying samples
            temp.rename(columns={"Name":"cellID"}, inplace=True)
            
            pos_cols = [col for col in temp.columns if "Positivity" in col and "(MV" in col]            # identify columns that says if protein marker is present (1) or abscent (0)
            data = temp[~(temp[pos_cols].eq(-1).any(axis=1))]                                           # filter out those cells that does contain a NaN value (-1) in the previous columns
            columns_nuc = data.columns[                                                                 # keep columns with protein intensity in the nucleus
                (data.columns.isin(['cellID']))|
                ((data.columns.str.contains("MV - NUC - "))&(~data.columns.str.contains("Type")))
            ]
            columns_pos = data.columns[data.columns.isin(['cellID', 'X-coordinate', 'Y-coordinate'])]   # keep spatial columns
            
            adata = AnnData(                                                                            # generate AnnData file to include spatial data
                data[columns_nuc].iloc[:,1:],
                obsm={
                    "spatial": data[columns_pos].iloc[:,1:].to_numpy(),
                    "ID_cell":data[['cellID']].to_numpy(),
                    "Positivity":data[pos_cols].to_numpy()
                }
            )
            adata_dicts[f"{data.iloc[0,0][:2]}"] = adata
            data_dicts[f"{data.iloc[0,0][:2]}"] = data
        
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
        if all(row[m] == v for m, v in rule.items()):
            return str(cell_type)
    return "Other cells"

## Labeling cell types
def labeling_cell_types(data_dicts, adata_dicts, cell_type_dict, output_dir, save_anndata=True):
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
    output_dir : Path
        Path to the output directory.
    save_anndata : Bool (Optional; Default is True)
        Boolean value to decide if AnnData files will be saved in anndata directory or not.

    Return
    ------
    adata_dicts : Dict[str:AnnData]
       Dictionary containing multiple updated annotated data matrices. 
    """
    for (k_data, data), (k, adata) in zip(data_dicts.items(), adata_dicts.items()):
        adata.obs['clusters'] = data.apply(assign_cell_type, axis=1, args=(cell_type_dict,)).tolist()
        for cell_type, _ in cell_type_dict.items():
            adata.obs[f"only {cell_type}"] = [t if t==cell_type else "Other cells" for t in adata.obs['clusters']]
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
    adata_dicts = load_anndata_files(output_dir)
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
        adata_dicts = labeling_cell_types(data_dicts, adata_dicts, cell_type_dict, output_dir, save_anndata = config['locally_save_anndata_files'])
        if config['locally_save_anndata_files']: print("Anndata saved!")
        else: print("Anndata not saved...")
    else: print("Anndata loaded!")
    return adata_dicts


## Plotting spacial data
def plot_spatial(adata_dicts,custom_colors,output_dir,overwrite_existing_files=False,dpi=300,size=100):
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
    size : int (Optional; Default is 100)
        Scatter dots size for spatial plot.
    """
    plots_path = output_dir / "plots" 
    plots_path.mkdir(parents=True, exist_ok=True)
    for k, adata in adata_dicts.items():
        title_name = f"Sample {k} ({adata.n_obs} cells)"
        save_namefile = plots_path / f"Spatial - {title_name}.png"
        if save_namefile.exists() and not overwrite_existing_files:
            print(f"File 'Spatial - {title_name}.png' already exists...")
            for cell_type, selected_color in custom_colors.items():
                if cell_type=='Other cells': continue
                if (adata.obs['clusters'] == cell_type).sum()==0: continue
                
                title_name2 = f"Sample {k} - {cell_type} ({adata.n_obs} cells)"
                save_namefile2 = plots_path / f"Spatial - {title_name2}.png"
                if save_namefile2.exists() and not overwrite_existing_files:
                    print(f"File 'Spatial - {title_name2}.png' already exists...")
                    continue
                color_spatial = f"only {cell_type}"
                cmap_gene = selected_color
                own_palette_list = [selected_color,custom_colors['Other cells']] if cell_type<'Other cells' else [custom_colors['Other cells'],selected_color]
                own_palette = ListedColormap(own_palette_list)
                cmap_gene = ListedColormap(cmap_gene)
                
                fig, ax = plt.subplots()
                sq.pl.spatial_scatter(
                    adata, 
                    shape=None, 
                    color=color_spatial,
                    title=title_name2,
                    dpi=dpi,
                    cmap=cmap_gene,
                    palette=own_palette,
                    size=size,
                    ax=ax)
                ax.set_facecolor("black")
                ax.set_aspect('equal')
                fig.tight_layout()
                plt.savefig(save_namefile2)
                plt.close()
                print(f"File 'Spatial - {title_name2}.png' created!")
            continue
        color_spatial = "clusters"
        cmap_gene = None
        which_colors = []
        for cell_type,selected_color in custom_colors.items():
            if (adata.obs['clusters'] == cell_type).sum()!=0: which_colors.append(selected_color)
        own_palette_list = (which_colors)
        own_palette = ListedColormap(own_palette_list)
        
        fig, ax = plt.subplots()
        sq.pl.spatial_scatter(
            adata, 
            shape=None, 
            color=color_spatial,
            title=title_name,
            dpi=dpi,
            cmap=cmap_gene,
            palette=own_palette,
            size=size,
            ax=ax)
        ax.set_facecolor("black")
        ax.set_aspect('equal')
        fig.tight_layout()
        plt.savefig(save_namefile)
        plt.close()
        print(f"File 'Spatial - {title_name}.png' created!")
        
        for cell_type, selected_color in custom_colors.items():
            if cell_type=='Other cells': continue
            if (adata.obs['clusters'] == cell_type).sum()==0: continue
            
            title_name2 = f"Sample {k} - {cell_type} ({adata.n_obs} cells)"
            save_namefile2 = plots_path / f"Spatial - {title_name2}.png"
            if save_namefile2.exists() and not overwrite_existing_files: 
                print(f"File 'Spatial - {title_name2}.png' already exists...")
                continue
            color_spatial = f"only {cell_type}"
            cmap_gene = selected_color
            own_palette_list = [selected_color,custom_colors['Other cells']] if cell_type<'Other cells' else [custom_colors['Other cells'],selected_color]
            own_palette = ListedColormap(own_palette_list)
            cmap_gene = ListedColormap(cmap_gene)
            
            fig, ax = plt.subplots()
            sq.pl.spatial_scatter(
                adata, 
                shape=None, 
                color=color_spatial,
                title=title_name2,
                dpi=dpi,
                cmap=cmap_gene,
                palette=own_palette,
                size=size,
                ax=ax)
            ax.set_facecolor("black")
            ax.set_aspect('equal')
            fig.tight_layout()
            plt.savefig(save_namefile2)
            plt.close()
            print(f"File 'Spatial - {title_name2}.png' created!")

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
        title_name = f"Sample {k} ({adata.n_obs} cells)"
        save_namefile = output_dir / "results" / f"Cell type proportions - {title_name}.csv"
        if save_namefile.exists() and not overwrite_existing_files: 
            print(f"File 'Cell type proportions - {title_name}' already exists...")
            continue
        sum_dict = {}
        sum_dict["Total cells"] = adata.n_obs
        if sum_dict["Total cells"]==0: 
            print(f"No cells in '{title_name}' file. Skipping...")
            continue
        prop_dict = {}
        for cell_type,_ in custom_colors.items():
            sum_dict[cell_type]  = (adata.obs['clusters'] == cell_type).sum()
            prop_dict[cell_type] = str(np.round(sum_dict[cell_type]/sum_dict["Total cells"]*100,2))+"%"
        
        df = pd.DataFrame(data = {cell_type:[sum_dict[cell_type],prop_dict[cell_type]] for cell_type,_ in custom_colors.items()}, index = ["Total cells in cell type", "Percentage"])
        df.to_csv(save_namefile, index=False)
        print(f"File 'Cell type proportions - {title_name}.csv' created!")









