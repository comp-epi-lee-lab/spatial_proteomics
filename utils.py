
# Load needed libraries

from pathlib import Path
import os

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

## Loading adata files
def load_anndata_files():
    sample_names = pd.read_csv(f"{results_dir}Samples Id.csv")['Samples Id'].tolist()
    adata_dicts = {}
    for sample in sample_names:
        adata_dicts[sample] = sc.read_h5ad(f"{adata_dir}Sample_{sample}_after_umap.h5ad")
    return adata_dicts

## Obtaining filenames
def filenames():
    path = Path(files_dir)
    entries = list(path.iterdir())
    filenames = [entry for entry in entries if entry.is_file() and entry.match(f"*objects.{filetype}")]
    return filenames

## Cleaning data
def cleaned_data(file_names, results_dir, filetype):
    adata_dicts = {}
    data_dicts  = {}
    for filename in file_names:
        separator = '\t' if filetype=='tsv' else ','
        data = pd.read_csv(f"{filename}",sep=separator)
        data = data[data['Positivity - DAPI (MV - NUC)']==1]                                        # filter out cells without nucleus
        
        temp = pd.concat([data.iloc[:,3], data.iloc[:, 12:]], axis=1)                               # retain only "Name" and data columns
        temp['Name'] = [f"{s[-5]}{s[-2]}_{i:05}" for s,i in zip(temp['Name'],temp['Name'].index)]   # retain letter and number for identifying samples
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
    df.to_csv(f"{results_dir}Samples Id.csv", index=False)
    df = pd.DataFrame({"Positivity column names": list(data[pos_cols].columns)})
    df.to_csv(f"{results_dir}Positivity column names.csv", index=False)
    return data_dicts, adata_dicts

## Assign cell types
def assign_cell_type(row):
    for cell_type, rule in cell_types.items():
        if all(row[m] == v for m, v in rule.items()):
            return str(cell_type)
    return "Other cells"

## Labeling cell types
def labeling_cell_types(data_dicts, adata_dicts, cell_type_dict):
    for (k_data, data), (k, adata) in zip(data_dicts.items(), adata_dicts.items()):
        adata.obs['clusters'] = data.apply(assign_cell_type, axis=1).tolist()
        for cell_type, _ in cell_type_dict.items():
            adata.obs[f"only {cell_type}"] = [t if t==cell_type else "Other cells" for t in adata.obs['clusters']]
        adata_dicts[k] = adata
    return adata_dicts

## Saving adata files
def save_anndata_files(adata_dicts):
    for k, adata in adata_dicts.items():
        adata.write(f"{adata_dir}Sample_{k}_after_umap.h5ad")

## Plotting spacial data
def plot_spatial(adata_dicts,custom_colors,dpi=300,size=100):
    for k, adata in adata_dicts.items():
        title_name = f"Sample {k} ({adata.n_obs} cells)"
        save_namefile = f"{plots_dir}Spatial - {title_name}.png"
        if os.path.exists(save_namefile) and not override_existing_plots:
            for cell_type, selected_color in custom_colors.items():
                if cell_type=='Other cells': continue
                if (adata.obs['clusters'] == cell_type).sum()==0: continue
                
                title_name = f"Sample {k} - {cell_type} ({adata.n_obs} cells)"
                save_namefile = f"{plots_dir}Spatial - {title_name}.png"
                if os.path.exists(save_namefile) and not override_existing_plots: continue
                
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
        
        for cell_type, selected_color in custom_colors.items():
            if cell_type=='Other cells': continue
            if (adata.obs['clusters'] == cell_type).sum()==0: continue
            
            title_name = f"Sample {k} - {cell_type} ({adata.n_obs} cells)"
            save_namefile = f"{plots_dir}Spatial - {title_name}.png"
            if os.path.exists(save_namefile) and not override_existing_plots: continue
            
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

## Calculating cell proportions and saving them in csv files
def calculate_cell_proportions(adata_dicts,custom_colors,pie_plot=False):
    for k, adata in adata_dicts.items():
        title_name = f"Sample {k} ({adata.n_obs} cells)"
        save_namefile = f"{results_dir}Cell type proportions - {title_name}.csv"
        if os.path.exists(save_namefile) and not override_existing_plots: continue
            
        sum_dict = {}
        sum_dict["Total cells"] = adata.n_obs
        if sum_dict["Total cells"]==0: continue
        prop_dict = {}
        for cell_type,_ in custom_colors.items():
            sum_dict[cell_type]  = (adata.obs['clusters'] == cell_type).sum()
            prop_dict[cell_type] = str(np.round(sum_dict[cell_type]/sum_dict["Total cells"]*100,2))+"%"
        
        df = pd.DataFrame(data = {cell_type:[sum_dict[cell_type],prop_dict[cell_type]] for cell_type,_ in custom_colors.items()}, index = ["Total cells in cell type", "Percentage"])
        df.to_csv(f"{results_dir}Cell type proportions - {title_name}.csv", index=False)
        
        if pie_plot: plot_pie_plots(sum_dict,custom_colors,title_name)