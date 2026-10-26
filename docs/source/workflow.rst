Analysis workflow
=================

Common preprocessing
--------------------

All analysis modes begin by:

#. reading Visiopharm object-level tables;
#. constructing per-sample AnnData objects;
#. applying user-defined rules for major cell types; and
#. storing cell-level annotations and spatial/intensity information.

Execution then branches according to ``type_of_analysis``.

Supervised branch
-----------------

The supervised branch applies user-defined subtype rules, calculates
cell-population counts and proportions, generates hierarchy information, and
produces requested spatial plots.

Unsupervised branch
-------------------

The unsupervised branch:

#. creates cell-type proportion matrices;
#. identifies binary protein marker patterns observed in the data;
#. quantifies those patterns within selected major cell types;
#. stores discovered phenotype labels/patterns in AnnData;
#. optionally compares phenotypes between two sample groups; and
#. optionally summarizes protein marker intensities within discovered phenotypes.
#. optionally summarizes major cell types neighborhood analysis.
#. optionally describes positive markers within discovered phenotypes using sample-relative low/medium/high intensity states.

Full branch
-----------

Full mode runs the unsupervised branch first, converts discovered protein marker
patterns into subtype definitions, and then invokes the supervised workflow.

Questions answered by each mode
-------------------------------

**Supervised**
   Which predefined phenotypes are present, how abundant are they, and how are
   they spatially distributed?

**Unsupervised**
   Which marker combinations are observed within a selected parent population,
   how abundant are they, and do their proportions differ between groups?

   Where is protein abundance located? Which major cell populations spatially associate?
   How strong is protein marker expression inside an already-discovered phenotype?

Neither mode alone establishes biological lineage or causality.
