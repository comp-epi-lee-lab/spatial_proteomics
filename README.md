# Spatial Proteomics Analysis (SpPrAn)

**SpPrAn (Spatial Proteomics Analysis)** is a Python-based pipeline for analyzing single-cell spatial proteomics data exported from Visiopharm®. It classifies cells using user-defined protein-marker positivity rules, supports hierarchical cell subtypes, generates spatial phenotype maps, and reports cell-population counts and proportions.

SpPrAn is intended for research-grade analyses in cancer biology, tumor-microenvironment studies, and other academic or translational spatial proteomics applications.

## Key features

* Process multiple Visiopharm® object-level `.tsv` or `.csv` files in a single run.
* Define major cell populations using positive, negative, or unconstrained protein markers.
* Create hierarchical cell subtypes using additional marker rules.
* Store annotated single-cell data as AnnData (`.h5ad`) files for downstream analysis.
* Generate overview, highlighted, or custom spatial phenotype plots.
* Calculate cell counts and percentages relative to both parent populations and all cells.
* Visualize marker-defined cell-type and subtype hierarchies.
* Create or edit YAML configuration files using the browser-based **SpPrAn Configuration Builder**.
* Run the complete pipeline either from the command line or through the desktop graphical interface.

## Documentation

Full documentation is maintained in the `docs/` directory and is intended for publication through **Read the Docs**.

The documentation includes installation, command-line and desktop-app usage, configuration creation, use of the SpPrAn Configuration Builder, marker-based cell-type/subtype definitions, output interpretation, and Python API documentation.

Documentation preview version is available [here](http://htmlpreview.github.io/?https://github.com/comp-epi-lee-lab/spatial_proteomics/blob/main/docs/_build/html/index.html).

## Ways to use SpPrAn

### 1. Desktop application

Pre-built desktop applications for **Windows** and **macOS** are available through (GitHub Releases)[https://github.com/comp-epi-lee-lab/spatial_proteomics/releases]. The current software releases are maintained by [Dr. Sergio Zamora-Erazo](mailto:sergio.zamorae@gmail.com).

The desktop application is intended for users who prefer not to install Python or work from the command line. It lets users create, edit, or open a configuration file, select input/output directories, launch the pipeline, and view runtime messages.

The desktop application includes the **SpPrAn Configuration Builder**, which opens locally in the user's default web browser.

> The Configuration Builder operates locally in the browser. Representative data files and folders selected through the builder are used only to read local column names or filenames needed for configuration; they are not uploaded to an external server.

#### Download and startup

Choose the build that matches your operating system and processor architecture:

* **Windows — x86_64**
  For 64-bit Intel or AMD processors.

* **Windows — ARM64**
  For ARM-based Windows devices, including Copilot+ PCs.

* **macOS — Apple Silicon**
  For Macs using M1, M2, M3, or M4 processors.

### Windows

1. Download the appropriate Windows ZIP file from the latest GitHub Release.
2. Extract the contents of the ZIP file.
3. Open the extracted folder.
4. Double-click the SpPrAn GUI `.exe` file to launch the application.

Windows may display a SmartScreen or security warning if the application is not digitally signed. Before continuing, verify that the application was downloaded from this project's official GitHub Releases page.

### macOS

1. Download the macOS Apple Silicon ZIP file from the latest GitHub Release.
2. Extract the ZIP file.
3. Open the extracted folder.
4. Double-click the SpPrAn GUI `.app` file to launch the application.

If macOS prevents the application from opening because the developer cannot be verified, open **System Settings → Privacy & Security** and choose **Open Anyway** for SpPrAn GUI.

Only bypass this warning if the application was downloaded directly from this project's official GitHub Releases page.

### 2. Python / command line

Users who prefer direct access to the Python environment can clone and install the repository.

## Requirements

* Python >= 3.12

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

The builder can create a new configuration, restore/edit an existing YAML file, identify positivity columns from a representative `.tsv`/`.csv`, define major populations and nested subtypes, assign colors, configure AnnData/output behavior, configure spatial plotting, and preview/download the YAML file.

### Manual YAML configuration

A current template is provided at:

```text
config/config_example.yaml
```

Create a working copy:

```bash
cp config/config_example.yaml config/config.yaml
```

Then edit the copy for the intended dataset.

## Running the pipeline

```bash
sppran --config config/config.yaml
```

Outputs are written to the directory specified under `workspace.output_dir`.

## Input data

SpPrAn is designed for object-level spatial proteomics tables exported from Visiopharm®. Input files are expected to contain one row per segmented cell/object and the positivity columns referenced by the configuration file.

The biological interpretation of a SpPrAn population depends on the marker panel, staining performance, segmentation strategy, and positivity thresholds used upstream.

## Expected outputs

Depending on the configuration, SpPrAn can generate:

* annotated AnnData (`.h5ad`) files;
* primary cell-type annotations;
* hierarchical subtype annotations;
* cell-count and population-proportion tables;
* percentages relative to parent populations and to all cells;
* marker-definition and hierarchy information;
* spatial overview plots;
* population-highlighted spatial plots; and
* cell-type/subtype hierarchy visualizations.

Detailed output descriptions and interpretation guidance are provided in the documentation.

## Interpretation

SpPrAn assigns **marker-defined phenotypes**. These labels should not be interpreted as independent biological ground truth.

For example:

```text
Immune cells
└── T cells
    └── CD8 T cells
```

represents nested marker-based phenotype definitions, not by itself a developmental lineage or differentiation trajectory.

## Citation

A formal software citation will be provided with the associated SpPrAn publication.

If SpPrAn is used in academic work before that citation is finalized, please contact the project authors ([Dr. Janusz Franco-Barraza](mailto:Janusz.FrancoBarraza@fccc.edu) for biological matters, [Dr. Hayan Lee](mailto:Hayan.Lee@fccc.edu) for computational matters, and [Dr. Sergio Zamora-Erazo](mailto:sergio.zamorae@gmail.com) for technical matters) for the preferred citation.

## License

**License pending**

No open-source license has currently been granted for SpPrAn. Until licensing terms are formally specified, the source code remains subject to applicable copyright restrictions.

The repository is publicly available for scientific transparency, evaluation, and reproducibility. Please contact the repository owner ([Dr. Hayan Lee](mailto:Hayan.Lee@fccc.edu)) before reuse, modification, redistribution, or incorporation into another software product.

## Contributing and issues

Bug reports, questions, and suggestions can be submitted through the GitHub Issues page.

## Contact

**Sergio Zamora-Erazo, Ph.D.**
<!-- Cukierman Lab, Greenberg Pancreatic Cancer Institute -->
<!-- Lee Lab, Cancer Epigenetics Institute -->
<!-- Fox Chase Cancer Center -->
Philadelphia, PA, USA

[sergio.zamorae@gmail.com](mailto:sergio.zamorae@gmail.com)
