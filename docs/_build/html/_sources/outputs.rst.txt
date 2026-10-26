Outputs
=======

Outputs depend on the analysis mode.

Shared outputs
--------------

AnnData objects can retain major cell-type labels together with spatial and
intensity data.

Supervised outputs
------------------

These can include:

* primary and subtype annotations;
* cell-count/proportion tables;
* percentages relative to parent populations and all cells;
* phenotype-definition tables;
* hierarchy visualizations; and
* spatial overview, highlight, or custom plots.

Unsupervised outputs
--------------------

Cell-type proportion matrices
   Sample-by-population tables summarizing major cell-type proportions.

Pie plots
   Optional sample-level visual summaries of population composition.

Discovered phenotype tables
   Protein-marker-pattern phenotypes observed within selected parent populations,
   together with within-parent cell counts and proportions.

Phenotype-proportion matrices
   Sample-by-phenotype tables for each explored parent.

AnnData phenotype annotations
   Discovered phenotype host names are stored in ``adata.obs`` and pattern
   definitions in ``adata.uns``.

Two-group comparison
   When two groups are configured, SpPrAn reports common and group-specific
   phenotypes and descriptive abundance summaries.

Statistical results
   Common phenotype proportions can be compared with two-sided Mann--Whitney U
   and Welch t-tests. Benjamini--Hochberg FDR-adjusted p-values are exported.

Intensity outputs
   For selected markers, SpPrAn exports per-phenotype Tukey box-plot summaries,
   complete intensity-data tables, and protein-intensity spatial maps.

(**New in v0.3.1**) Neighborhood analysis outputs
   When selected, SpPrAn exports neighborhood-enrichment and interaction-matrix 
   overview heatmaps, and distance-dependent co-occurrence plots per sample.

(**New in v0.3.1**) Advanced phenotype classification
   When selected, SpPrAn reports:
   
   - Low/medium/high intensity abundance for positive markers within each
   phenotype.
   - Observed combinations of marker-intensity states and their frequencies.
   - Phenotype-level summary of intensity-profile complexity.

Full-mode outputs
-----------------

Full mode produces exploratory outputs and then adds hierarchy-oriented
supervised outputs based on the discovered phenotype definitions.

Interpretation
--------------

Phenotype proportions are conditional on the parent population used for
discovery. For example, a phenotype representing 20% of T cells does not necessarily
represent 20% of all cells in the tumor.

Neighborhood outputs describe spatial association rather than causal
communication.

Advanced intensity classes are relative within each sample and should not be
interpreted as absolute cross-sample expression categories.

Statistical results should be interpreted together with detection frequency,
effect magnitude, biological sample size, multiple-testing adjustment, and
study design.
