SpPrAn documentation
====================

**SpPrAn (Spatial Proteomics Analysis)** is a Python pipeline for analyzing
single-cell spatial proteomics data exported from Visiopharm®. Cells are
classified using user-defined protein-marker rules and can be organized into
hierarchical cell populations and subtypes.

SpPrAn can be run from Python/the command line or from a packaged desktop
application for users who prefer a graphical workflow.

Version 0.3.0
-------------

SpPrAn 0.3.0 introduces separate supervised and unsupervised analysis
workflows, plus an optional full phenotype-discovery workflow that connects
exploratory phenotype discovery with the established supervised analysis.

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
