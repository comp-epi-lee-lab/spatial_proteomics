from .utils import *

class Config:
    def __init__(self, workspace_dir, files_dir, filetype):
        self.workspace_dir = workspace_dir
        self.files_dir = files_dir
        self.filetype = filetype
        self.results_dir = f"{workspace_dir}results/"
        self.plots_dir = f"{workspace_dir}plots/"
        self.adata_dir = f"{workspace_dir}adata/"

