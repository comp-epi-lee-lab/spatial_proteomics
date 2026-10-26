Running SpPrAn
==============

Command line
------------

.. code-block:: bash

   sppran --config config/config.yaml

Desktop application
-------------------

The graphical application provides the same underlying analysis pipeline with a
user interface for configuration selection, workspace selection, execution,
and runtime messages. See :doc:`desktop_app`.

Choosing a mode
---------------

Use **Supervised** for a hypothesis-driven analysis with predefined subtypes.

Use **Unsupervised** for exploratory marker-pattern discovery within known
major parent populations.

Enable **full phenotype discovery** when the discovered patterns should be
converted into subtype definitions and passed through the supervised workflow.

Typical supervised workflow
---------------------------

#. Export object-level measurements from Visiopharm®.
#. Create or update a YAML configuration.

   * Define biologically appropriate protein marker rules.
   * Add hierarchical subtypes where supported by the marker panel.
   * Select spatial plotting/output options.

#. Run SpPrAn from the desktop application or command line.
#. Review population counts, hierarchy reports, and spatial maps.

Typical unsupervised workflow
-----------------------------

#. Export object-level measurements from Visiopharm®.
#. Create or update a YAML configuration.

   * Define biologically appropriate protein marker rules.
   * Select populations to explore.
   * Select exploratory protein markers and their corresponding short names.
   * Optionally enable full phenotype discovery option.
   * Optionally define two sample groups for all the samples in the run.
   * Optionally select protein markers for intensity summaries.
   * Optionally enable neighborhood analysis.
   * Optionally enable advanced phenotype classification.

#. Run SpPrAn from the desktop application or command line.
#. Review phenotype patterns and within-parent proportions.
#. If groups were supplied, interpret descriptive differences, detection 
   frequency, nominal p-values, and FDR-adjusted results together.
#. If intensity protein markers were supplied, review boxplot summaries and 
   outliers from intensities by phenotype in each selected cell type.