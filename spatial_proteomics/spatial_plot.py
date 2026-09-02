from matplotlib import rcParams
import pandas as pd
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt

import warnings
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
from spatial_proteomics.utils import export_legend_png

import logging

logger = logging.getLogger(__name__)

def plotting_helper(custom_colors, overwrite_existing_files, selected_plots, dpi, size, ncol, plots_path, k, adata, folder_name):
    if folder_name == 'General':
        universe_to_color = 'clusters' 
        title_name = f"Sample {k} ({adata.n_obs} cells)"
    else:
        universe_to_color = f"Subtypes of {folder_name}"
        if folder_name in adata.obs['clusters'].unique():
            total_subcells = (adata.obs['clusters'] == folder_name).sum()
        else:
            total_subcells = 0
            for idx in adata.obs.columns[adata.obs.columns.str.contains("Subtypes of")].unique():
                if folder_name in adata.obs[idx].unique():
                    total_subcells = (adata.obs[idx] == folder_name).sum()
                    break
        if total_subcells == 0: 
            logger.warning(f"No cells were found for {folder_name} in Sample {k}. Skipping plotting for this cell type...")
            return
        title_name = f"Sample {k} - {folder_name} ({total_subcells} cells)"
    save_namefile = plots_path / f"Spatial - {title_name}.png"

    if save_namefile.exists() and not overwrite_existing_files:
        if selected_plots != "highlight":
            logger.warning(f"File 'Spatial - {title_name}.png' already exists...")
        if selected_plots != "overview":
            for cell_type, selected_color in custom_colors.items():
                if cell_type=='Other cells': continue
                if cell_type.startswith('Other'): continue
                if cell_type not in adata.obs[universe_to_color].values: continue
                if (adata.obs[universe_to_color] == cell_type).sum()==0: continue
                title_name2 = f"Sample {k} - {cell_type} ({adata.n_obs} cells)" if folder_name == 'General' else f"Sample {k} - {folder_name} - {cell_type} ({total_subcells} cells)"
                check_path = plots_path / f"Other {folder_name} plots"
                check_path.mkdir(parents=True, exist_ok=True)
                save_namefile2 = check_path /  f"Spatial - {title_name2}.png"
                if save_namefile2.exists() and not overwrite_existing_files:
                    logger.warning(f"File 'Spatial - {title_name2}.png' already exists...")
                    continue
                # adata.obs[f"only {cell_type}"]
                # adata.obs[f"Subtypes of {cell_type}"] # Sub universe
                color_spatial = f"Only {cell_type}" if cell_type == folder_name or folder_name == 'General' else f"only {cell_type}"
                color_to_check = f"{folder_name} - only {cell_type}"
                color_spatial = color_to_check if adata.obs.columns.isin([color_to_check]).any() else color_spatial
                    # cmap_gene = selected_color
                # own_palette_list = [selected_color,custom_colors['Other cells']] if cell_type<'Other cells' else [custom_colors['Other cells'],selected_color]
                own_palette_list = [custom_colors[label] for label in sorted([cell_type, f"Other cells"])] if folder_name=='General' or (folder_name!='General' and not any(adata.obs[color_spatial]==f"Other {folder_name}")) else [custom_colors[label] for label in sorted([cell_type, f"Other {folder_name}", f"Other cells"])]
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
                            # legend_loc="lower center",
                            ax=ax)
                ax.set_facecolor("black")
                ax.set_aspect('auto')
                handles, labels = ax.get_legend_handles_labels()
                check_path_legend2 = check_path / "Plot Legends"
                check_path_legend2.mkdir(parents=True, exist_ok=True)
                save_namefile_legend2 = check_path_legend2 /  f"Spatial - Legend - {title_name2}.png"
                export_legend_png(ax,save_namefile_legend2,handles, labels, horizontal=False)
                if ax.get_legend() is not None: ax.get_legend().remove()
                ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.03), ncol=ncol, frameon=False,)
                fig.subplots_adjust(bottom=0.1)
                pos = ax.get_position()
                ax.set_position([pos.x0,pos.y0+0.05,pos.width,pos.height])
                # ax.set_aspect('equal')
                # fig.subplots_adjust(right=0.8)
                # fig.tight_layout()
                plt.savefig(save_namefile2, dpi=dpi)
                plt.close(fig)
        return
    if selected_plots != "highlight":
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
                    # legend_loc="lower center",
                    ax=ax)
        ax.set_facecolor("black")
        ax.set_aspect('auto')
        handles, labels = ax.get_legend_handles_labels()
        check_path_legend = plots_path / "Plot Legends"
        check_path_legend.mkdir(parents=True, exist_ok=True)
        save_namefile_legend = check_path_legend /  f"Spatial - Legend - {title_name}.png"
        export_legend_png(ax,save_namefile_legend,handles, labels, horizontal=False)
        if ax.get_legend() is not None: ax.get_legend().remove()
        ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.03), ncol=ncol, frameon=False,)
        fig.subplots_adjust(bottom=0.1)
        pos = ax.get_position()
        ax.set_position([pos.x0,pos.y0+0.05,pos.width,pos.height])
        # ax.set_aspect('equal')
        # fig.subplots_adjust(right=0.8)
        # plt.legend(ncol=1,loc='center right',bbox_to_anchor=(0.5, -0.05))
        # fig.tight_layout()
        plt.savefig(save_namefile, dpi=dpi)
        plt.close(fig)
        # print(f"File 'Spatial - {title_name}.png' created!")
        # print(f"Spatial plots for the cell types in {folder_name} has been created!")
        
    if selected_plots != "overview":
        for cell_type, selected_color in custom_colors.items():
            if cell_type=='Other cells': continue
            if cell_type.startswith('Other'): continue
            if cell_type not in adata.obs[universe_to_color].values: continue
            if (adata.obs[universe_to_color] == cell_type).sum()==0: continue
                
            title_name2 = f"Sample {k} - {cell_type} ({adata.n_obs} cells)" if folder_name == 'General' else f"Sample {k} - {folder_name} - {cell_type} ({total_subcells} cells)"
            check_path = plots_path / f"Other {folder_name} plots"
            check_path.mkdir(parents=True, exist_ok=True)
            save_namefile2 = check_path / f"Spatial - {title_name2}.png"
            if save_namefile2.exists() and not overwrite_existing_files: 
                logger.warning(f"File 'Spatial - {title_name2}.png' already exists...")
                continue
            color_spatial = f"Only {cell_type}" if cell_type == folder_name or folder_name == 'General' else f"only {cell_type}"
            color_to_check = f"{folder_name} - only {cell_type}"
            color_spatial = color_to_check if adata.obs.columns.isin([color_to_check]).any() else color_spatial
                # cmap_gene = selected_color
            # own_palette_list = [selected_color,custom_colors['Other cells']] if cell_type<'Other cells' else [custom_colors['Other cells'],selected_color]
            own_palette_list = [custom_colors[label] for label in sorted([cell_type, f"Other cells"])] if folder_name=='General' or (folder_name!='General' and not any(adata.obs[color_spatial]==f"Other {folder_name}")) else [custom_colors[label] for label in sorted([cell_type, f"Other {folder_name}", f"Other cells"])]
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
                        # legend_loc="lower center",
                        ax=ax)
            ax.set_facecolor("black")
            ax.set_aspect('auto')
            handles, labels = ax.get_legend_handles_labels()
            check_path_legend2 = check_path / "Plot Legends"
            check_path_legend2.mkdir(parents=True, exist_ok=True)
            save_namefile_legend2 = check_path_legend2 /  f"Spatial - Legend - {title_name2}.png"
            export_legend_png(ax,save_namefile_legend2,handles, labels, horizontal=False)
            if ax.get_legend() is not None: ax.get_legend().remove()
            ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.03), ncol=ncol, frameon=False,)
            fig.subplots_adjust(bottom=0.1)
            pos = ax.get_position()
            ax.set_position([pos.x0,pos.y0+0.05,pos.width,pos.height])
            # ax.set_aspect('equal')
            # fig.subplots_adjust(right=0.8)
            # fig.tight_layout()
            plt.savefig(save_namefile2, dpi=dpi)
            plt.close(fig)

def plot_spatial(adata_dicts,custom_colors,output_dir,plot_options,custom_plot,overwrite_existing_files=False,dpi=300,size=50,ncol=4):
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
    samples_filepath = output_dir / "results" / "Samples Id.csv"
    if samples_filepath.exists(): sample_names_dict = pd.read_csv(samples_filepath,dtype={"Samples Id": str}).set_index('Samples Id').to_dict()['Samples filename']
    custom_plot_dict = {v['fileName'] : v['type'] for w,v in enumerate(custom_plot)}
    for k, adata in adata_dicts.items():
        if "color_keys" in adata.uns and "color_vals" in adata.uns:
            for key,val in zip(adata.uns['color_keys'], adata.uns['color_vals']):
                if key not in custom_colors.keys(): custom_colors[key] = val
        selected_plots = plot_options
        if plot_options == 'custom':
            if custom_plot_dict[sample_names_dict[k]] == 'none':
                continue
            selected_plots = custom_plot_dict[sample_names_dict[k]]
        for folder_name in adata.uns['folder_names']:
            plots_path_folder = plots_path / folder_name
            plots_path_folder.mkdir(parents=True, exist_ok=True)
            plotting_helper(custom_colors, overwrite_existing_files, selected_plots, dpi, size, ncol, plots_path_folder, k, adata, folder_name)
        logger.warning(f"All spatial plots for cell subtypes in Sample {k} have been created!")