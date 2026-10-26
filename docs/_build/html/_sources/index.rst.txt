SpPrAn documentation
====================

**SpPrAn (Spatial Proteomics Analysis)** is a Python pipeline for analyzing
single-cell spatial proteomics data exported from Visiopharm®. Cells are
classified using user-defined protein-marker rules and can be organized into
hierarchical cell populations and subtypes.

SpPrAn can be run from Python/the command line or from a packaged desktop
application for users who prefer a graphical workflow.

Version 0.3.1
-------------

SpPrAn 0.3.1 expands the unsupervised/full workflow with continuous
protein-intensity spatial maps, fixed-radius neighborhood analysis, and
intensity-aware advanced phenotype classification.

.. important::

   SpPrAn phenotypes are protein-marker-defined operational annotations. Their
   biological meaning depends on the antibody panel, tissue context,
   segmentation strategy, staining performance, and positivity thresholds.

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   installation
   desktop_app
   analysis_modes
   optional_analyses
   config_builder
   workflow
   configuration
   usage
   outputs
   release_notes

.. toctree::
   :maxdepth: 2
   :caption: Developer reference

   api
