Analysis modes
==============

SpPrAn provides two user-facing choices in the Configuration Builder:
**Supervised** and **Unsupervised**. A third configuration value, ``full``, is
used when optional full phenotype discovery is enabled.

Supervised mode
---------------

``type_of_analysis: supervised``

Supervised mode is hypothesis-driven. Users define major cell types and,
optionally, nested subtypes using protein-marker positivity rules.

Use this mode when the protein marker combinations defining the populations of interest
are known in advance.

The workflow includes:

* primary cell-type annotation;
* hierarchical subtype annotation;
* cell counts and proportions;
* phenotype-hierarchy outputs; and
* spatial overview, highlight, or custom plots.

Unsupervised mode
-----------------

``type_of_analysis: unsupervised``

Unsupervised mode keeps major cell types user-defined, but discovers
subphenotypes from observed combinations of selected exploratory protein markers within
one or more of those parent populations.

The term **unsupervised** does not mean that all cells are clustered without
prior biological information. Parent populations remain protein-marker-rule-defined.
The exploratory step enumerates observed binary protein marker patterns within them.

The workflow can generate:

* sample-by-cell-type proportion matrices;
* optional sample-level pie plots;
* observed phenotype labels and protein-marker-pattern mappings;
* phenotype proportions within selected parent populations;
* optional two-group phenotype comparisons; and
* optional protein-marker intensity summaries.

Version 0.3.1 additionally provides optional:

* continuous protein-marker intensity spatial maps;
* fixed-radius neighborhood analysis; and
* advanced low/medium/high phenotype-intensity characterization.

See :doc:`optional_analyses`.

Two-group comparison
~~~~~~~~~~~~~~~~~~~~

Samples can optionally be assigned to two groups, such as Control and
Treatment. For phenotypes found in both groups, SpPrAn can report descriptive
summaries and perform:

* two-sided Mann--Whitney U tests;
* Welch independent-samples t-tests; and
* Benjamini--Hochberg FDR correction.

Group-specific phenotypes are also reported.

.. caution::

   Statistical significance should be interpreted together with biological
   sample size, phenotype-detection frequency, effect magnitude, multiple-test
   correction, and experimental design.

Protein-marker-intensity summaries
~~~~~~~~~~~~~~~~~~~~~~~~~~

Selected protein-marker intensities can be summarized within discovered
phenotypes. SpPrAn exports Tukey-style box-plot summary values and complete intensity-data tables for
downstream visualization.

Full phenotype discovery
------------------------

``type_of_analysis: full``

The Configuration Builder exposes full phenotype discovery as an option in the
unsupervised workflow.

Full mode:

#. runs unsupervised protein-marker-pattern discovery;
#. combines the observed phenotype definitions across samples;
#. converts those protein marker patterns into SpPrAn subtype definitions;
#. regenerates annotated AnnData objects; and
#. continues through the supervised workflow.

.. caution::

   If the total of discovered phenotypes for a cell type is higher than 512, 
   they will be omitted from the supervised analysis.

Biological interpretation
-------------------------

A discovered phenotype is an observed combination of protein marker positivity states
within a defined parent population. It is not automatically a distinct cell
lineage, differentiation state, or functional state.

Interpret the phenotype using the tissue context and the known specificity and
function of the protein markers in the panel.
