# Spatial Proteomics Analysis (SpPrAn)

**SpPrAn (Spatial Proteomics Analysis)** is a Python-based pipeline for analyzing single-cell spatial proteomics data exported from Visiopharm®. It classifies cells using protein-marker positivity rules, supports hierarchical cell subtypes, provides exploratory phenotype discovery, generates spatial and statistical outputs, and reports cell-population counts and proportions.

SpPrAn is intended for research-grade analyses in cancer biology, tumor-microenvironment studies, and other academic or translational spatial proteomics applications.

**Current source version: 0.3.0**

## Key features

- Process multiple Visiopharm® object-level `.tsv` or `.csv` files in a single run.
- Define major cell populations using positive, negative, or unconstrained protein markers.
- Choose between **supervised** and **unsupervised** analysis workflows.
- In supervised mode, define hierarchical cell subtypes directly from protein-marker rules.
- In unsupervised mode, discover observed protein-marker-defined phenotypes within selected major cell types.
- Optionally compare discovered phenotypes between two sample groups.
- Optionally export protein-marker intensity summaries for discovered phenotypes.
- Optionally run **full phenotype discovery**, which converts discovered phenotypes into subtype definitions and continues through the supervised workflow.
- Store annotated single-cell spatial proteomics data as AnnData (`.h5ad`) files for downstream analysis.
- Generate overview, highlighted, custom spatial phenotype plots, and optional pie plots.
- Calculate cell counts and proportions.
- Visualize protein-marker-defined cell-type and subtype hierarchies.
- Create or edit YAML configuration files using the browser-based **SpPrAn Configuration Builder**.
- Run the pipeline either from the command line or through the desktop graphical interface.

## Analysis modes

SpPrAn 0.3.0 introduces separate analysis modes for hypothesis-driven and exploratory spatial-proteomics workflows.

### Supervised mode

In **supervised mode**, the user defines both major cell types and, when desired, nested subtypes using protein-marker presence/absence rules.

This mode is appropriate when the expected biological populations are already known from the marker panel or experimental hypothesis.

Supervised mode includes:

- protein-marker-rule-based cell annotation;
- hierarchical subtype assignment;
- spatial phenotype plotting;
- cell counts and proportions; and
- phenotype-hierarchy visualization.

### Unsupervised mode

In **unsupervised mode**, the user still defines the **major cell types** using protein-marker rules, but does not need to predefine all subtypes. Instead, SpPrAn explores observed combinations of selected exploratory protein markers within one or more selected major cell types.

This is best understood as **exploratory protein-marker-pattern discovery within biologically defined parent populations**. It is not unsupervised clustering in an embedding or latent feature space.

The unsupervised workflow can:

- generate cell-type proportion matrices;
- optionally create sample-level pie plots;
- identify observed binary protein-marker phenotypes within selected parent cell types;
- calculate phenotype proportions;
- optionally compare two user-defined sample groups;
- perform Mann–Whitney and Welch t-test comparisons with FDR-corrected results; and
- optionally export protein-marker-intensity summaries for downstream boxplot analysis.

### Full phenotype discovery

The Configuration Builder also provides an optional **full phenotype discovery** workflow. Internally, this uses `type_of_analysis: full`.

This workflow first performs unsupervised phenotype discovery, then converts the discovered protein-marker patterns into subtype definitions and continues with the standard supervised analysis workflow.

This option is useful when the goal is to discover protein-marker-defined phenotypes first and then obtain the downstream subtype-oriented outputs produced by supervised SpPrAn.

> **Biological interpretation:** discovered phenotypes are combinations of protein-marker positivity states observed within a user-selected parent population. They should be interpreted in the context of the antibody panel, tissue biology, segmentation strategy, and positivity thresholds. A discovered protein-marker combination is not automatically a distinct biological cell lineage or state.

## Documentation

SpPrAn documentation includes download and startup, command-line and desktop-app usage, configuration creation, analysis-mode selection, use of the SpPrAn Configuration Builder, protein-marker-based cell-type/subtype definitions, unsupervised phenotype discovery, output interpretation, and Python API documentation.

Documentation is available here:

https://sppran.readthedocs.io/en/latest/

## Ways to use SpPrAn

### 1. Desktop application

Pre-built desktop applications for **Windows** and **macOS** are available through [GitHub Releases](https://github.com/comp-epi-lee-lab/spatial_proteomics/releases). The software releases are maintained by [Dr. Sergio Zamora-Erazo](mailto:sergio.zamorae@gmail.com).

> The source code may be newer than the most recently published pre-built desktop release. Check the release version before downloading if exact version matching is important.

The desktop application is intended for users who prefer not to install Python or work from the command line. It lets users create, edit, or open a configuration file, select input/output directories, launch the pipeline, and view runtime messages.

The desktop application includes the **SpPrAn Configuration Builder**, which opens locally in the user's default web browser.

> The Configuration Builder operates locally in the browser. Representative data files and folders selected through the builder are used only to read local column names or filenames needed for configuration; they are not uploaded to an external server.

#### Download and startup

Choose the build that matches your operating system and processor architecture:

- **Windows — x86_64**  
  For 64-bit Intel or AMD processors.
- **Windows — ARM64**  
  For ARM-based Windows devices, including Copilot+ PCs.
- **macOS — Apple Silicon**  
  For Macs using M1, M2, M3, or M4 processors.

#### Windows

1. Download the appropriate Windows ZIP file from the latest GitHub Release.
2. Extract the contents of the ZIP file.
3. Open the extracted folder.
4. Double-click the SpPrAn GUI `.exe` file to launch the application.

Windows may display a SmartScreen or security warning if the application is not digitally signed. Before continuing, verify that the application was downloaded from this project's official GitHub Releases page.

#### macOS

1. Download the macOS Apple Silicon ZIP file from the latest GitHub Release.
2. Extract the ZIP file.
3. Open the extracted folder.
4. Double-click the SpPrAn GUI `.app` file to launch the application.

If macOS prevents the application from opening because the developer cannot be verified, open **System Settings → Privacy & Security** and choose **Open Anyway** for SpPrAn GUI.

Only bypass this warning if the application was downloaded directly from this project's official GitHub Releases page.

### 2. Python / command line

Users who prefer direct access to the Python environment can clone and install the repository.

## Requirements

- Python >= 3.12

Python dependencies are defined in `setup.py` and `requirements.txt`.

## Installation from source

```bash
git clone https://github.com/comp-epi-lee-lab/spatial_proteomics.git
cd spatial_proteomics
python -m pip install -e .
```

This installs the package and makes the `sppran` command-line entry point available.

## Configuration

SpPrAn analyses are controlled by a YAML configuration file.

### Configuration Builder

Open `config_builder.html` in a web browser, or launch it from the desktop SpPrAn application.

The builder can:

- create a new configuration or restore/edit an existing YAML file;
- select **Supervised** or **Unsupervised** mode;
- identify positivity columns from a representative `.tsv`/`.csv` file;
- define major cell populations;
- define hierarchical subtypes in supervised mode;
- select parent populations and exploratory protein markers in unsupervised mode;
- optionally configure two-group phenotype comparison;
- optionally select protein markers for intensity summaries;
- optionally enable full phenotype discovery;
- assign colors;
- configure AnnData/output behavior;
- configure spatial plotting where applicable; and
- preview/download the resulting YAML file.

### Manual YAML configuration

A current template is provided at:

```text
config/config_example.yaml
```

Create a working copy:

```bash
cp config/config_example.yaml selected/path/to/config.yaml
```

Then edit the copy for the intended dataset.

For backward compatibility, configurations without `type_of_analysis` are treated as supervised configurations.

## Running the pipeline

```bash
sppran --config path/to/config.yaml
```

Outputs are written to the directory specified under `workspace.output_dir`.

## Input data

SpPrAn is designed for object-level spatial proteomics tables exported from Visiopharm®. Input files are expected to contain one row per segmented cell/object and the positivity columns referenced by the configuration file.

The biological interpretation of a SpPrAn population depends on the marker panel, staining performance, segmentation strategy, and positivity thresholds used upstream.

## Expected outputs

Outputs depend on the selected analysis mode.

### Supervised outputs

These can include:

- annotated AnnData (`.h5ad`) files;
- primary cell-type annotations;
- hierarchical subtype annotations;
- cell-count and population-proportion tables;
- percentages relative to parent populations and to all cells;
- protein-marker-definition and hierarchy information;
- spatial overview plots;
- population-highlighted spatial plots; and
- cell-type/subtype hierarchy visualizations.

### Unsupervised outputs

These can include:

- annotated AnnData files containing discovered phenotype labels;
- cell-type proportion matrices;
- optional cell-type pie plots;
- phenotype-proportion tables for selected parent populations;
- protein-marker-pattern mappings for observed phenotypes;
- optional common/unique phenotype tables for two-group comparisons;
- optional Mann–Whitney and Welch t-test results with FDR correction;
- phenotype-proportion comparison plots; and
- optional protein-marker intensity summary tables for discovered phenotypes.

### Full phenotype-discovery outputs

Full mode first produces the exploratory outputs and then reuses discovered protein-marker patterns as subtype definitions for the supervised workflow.

Detailed output descriptions and interpretation guidance are provided in the documentation.

## Interpretation

SpPrAn assigns or discovers **protein-marker-defined phenotypes**. These labels should not be interpreted as independent biological ground truth.

For example:

```text
Immune cells
└── T cells
    └── CD8 T cells
```

represents nested protein-marker-based phenotype definitions, not by itself a developmental lineage or differentiation trajectory.

Likewise, a phenotype discovered in unsupervised mode represents an observed combination of positivity states among the selected exploratory protein markers within a defined parent population. Its biological meaning requires interpretation using the experimental protein-marker panel and cancer/tissue context.

## Citation

If you use **SpPrAn** in your research, please cite the software using the citation metadata provided in [`CITATION.cff`](CITATION.cff) after contacting the project authors ([Dr. Janusz Franco-Barraza](mailto:Janusz.FrancoBarraza@fccc.edu) and [Dr. Edna Cukierman](mailto:Edna.Cukierman@fccc.edu) for biological matters, and [Dr. Hayan Lee](mailto:Hayan.Lee@fccc.edu) and [Dr. Sergio Zamora-Erazo](mailto:sergio.zamorae@gmail.com) for bioinformatics affairs).

GitHub also provides formatted APA and BibTeX citations through the **Cite this repository** option on the repository page.

A dedicated SpPrAn methods/software publication is planned. Once that publication is available, it will be added as the preferred citation for the software.

### Publications using SpPrAn

SpPrAn has been used in published cancer research, including:

- Wang S, Lobello C, Gupta P, Francis R, Vistarop A, Shah B, **Zamora-Erazo S**, Ahmed F, Rushmore D, Kalgren A, Nelson D, Dupoux V, Basappa J, Baldwin DA, Fung H, Jillab M, **Lee H**, Cai KQ, Nejati R, Boldajipour B, **Franco-Barraza J**, Lee G, Nath K, Wasik MA, *Multipronged enhancement of CART cell activity against ROR1 neoantigen by inhibition of H3K27 methyltransferase EZH2.* **Cell Reports Medicine**. 2026. https://doi.org/10.1016/j.xcrm.2026.103060

In this publication, SpPrAn was used as part of the spatial single-cell proteomics analysis workflow.

Please note that publications listed here demonstrate applications of SpPrAn and are not necessarily the preferred citation for the software itself.

## License

**License pending**

No open-source license has currently been granted for SpPrAn. Until licensing terms are formally specified, the source code remains subject to applicable copyright restrictions.

The repository is publicly available for scientific transparency, evaluation, and reproducibility. Please contact [Dr. Hayan Lee](mailto:Hayan.Lee@fccc.edu) before reuse, modification, redistribution, or incorporation into another software product.

## Contributing and issues

Bug reports, questions, and suggestions can be submitted through the GitHub Issues page.

## Contact

**Sergio Zamora-Erazo, Ph.D.**  
Philadelphia, PA, USA

[sergio.zamorae@gmail.com](mailto:sergio.zamorae@gmail.com)
