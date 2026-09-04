SpPrAn documentation
====================

**SpPrAn (Spatial Proteomics Analysis)** is a Python pipeline for analyzing
single-cell spatial proteomics data exported from Visiopharm®. Cells are
classified using user-defined protein-marker rules and can be organized into
hierarchical cell populations and subtypes.

SpPrAn can be run from Python/the command line or from a packaged desktop
application for users who prefer a graphical workflow.

.. note::

   SpPrAn assigns marker-defined phenotypes. Biological specificity depends on
   the antibody panel, tissue context, segmentation strategy, staining
   performance, and upstream positivity thresholds.

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   installation
   desktop_app
   config_builder
   workflow
   configuration
   usage
   outputs

.. toctree::
   :maxdepth: 2
   :caption: Developer reference

   api
