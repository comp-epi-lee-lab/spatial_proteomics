# Spatial Proteomics Analysis (SpPrAn)
A bioinformatic pipeline to analyze spatial proteomics samples data obtained via Visiopharm® using cell types defined by presence or abscent of protein markers.

---
## Introduction
This repository contains an end-to-end python-based pipeline for **processing and analyzing presence or abscent of protein markers in single cell spatial proteomics data**. The pipeline is designed for reproducible analysis using multiple objects outputs as inputs. Supports different cancer types as well as configurable protein markers and cell types via a YAML configuration file.

The primary use case is research-grade data analysis in an academic or translational research setting.

---
## Requirements
- Python ≥ 3.10

Additional dependencies are listed in:
- `requirements.txt`

This pipeline has been tested on macOS 26.1 (M2, RAM 8 GB) and macOS 15.7.2 (M4, RAM 32 GB).

---
## Installation
At this stage, the pipeline is intended to be run directly from the repository.
```
git clone https://github.com/comp-epi-lee-lab/spatial_proteomics.git
cd spatial_proteomics
pip install -r requirements.txt
```
---

## Configuration
Pipeline parameters, input paths, and output locations are controlled through a YAML configuration file.
1. Create a working configuration file:
`cp config/config_example.yaml config/config.yaml`
2. Edit `config/config.yaml` to specify:
    - input data locations
    - output directory
    - analysis parameters

The file `config_example.yaml` serves as a documented template and should not be edited directly.

## Running the Pipeline

Run the full end-to-end pipeline using:
`python sppran --config config/config.yaml`
Outputs will be written to the directory specified in the configuration file.

---
## Expected Outputs
The pipeline may generate:
* AnnData files for easy generation of plots
* Cell type proportion tables
* Spatial plots featuring cell types

Detailed descriptions of outputs are provided in the documentation.

## Documentation
A step-by-step description of the workflow is available (or forthcoming) on protocols.io, including guidance on expected inputs, runtime, and interpretation of results.

---
## License
***License pending institutional review**

This software is currently under review by the authors’ institution to determine the appropriate licensing terms. Until a license is formally specified, all rights are reserved.

The code is shared for evaluation, review, and reproducibility of published research.

Please contact the corresponding author before any reuse, modification, or redistribution.

---
## Citation
If you use this pipeline for academic research, please cite:
Zamora-Erazo, S., Franco-Barraza, J., Lee, H. **SpPrAn**: a Bioinformatics pipeline for multiple Spatial Proteomics Analyses. Preprint, 2025.

(A formal software citation will be added once licensing is finalized.)

### Contact
For questions, issues, or collaboration inquiries, please contact:
- Sergio Zamora-Erazo     
sergio.zamora-erazo@fccc.edu     
Lee Lab, Cancer Epigenetics Institute     
Fox Chase Cancer Center     
Philadelphia, PA, USA     
