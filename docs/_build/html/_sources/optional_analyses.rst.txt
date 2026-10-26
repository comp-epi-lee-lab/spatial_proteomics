Optional spatial and phenotype analyses
=======================================

SpPrAn 0.3.1 adds three optional analyses to the unsupervised/full workflow.
They complement phenotype discovery rather than replacing it.

Protein-marker intensity spatial plots
--------------------------------------

The ``intensity_spatial_plot`` option creates a spatial map for each sample and
protein marker.

Each segmented cell remains at its tissue coordinate and is colored by the
continuous intensity stored in the AnnData matrix.

Configuration example:

.. code-block:: yaml

   intensity_spatial_plot:
     dpi: 150
     scatter_point_size: 50

Biological interpretation
~~~~~~~~~~~~~~~~~~~~~~~~~

These maps answer **where high or low measured protein intensity occurs in the
tissue**.

They are conceptually similar to displaying gene expression on a UMAP, except
that the x/y positions are true tissue coordinates rather than embedding
coordinates.

.. caution::

   The current implementation scales each sample/protein marker plot independently from
   its observed minimum and maximum intensity. Colors should therefore be
   interpreted within a panel and should not be treated as a quantitatively
   harmonized scale across samples.

Fixed-radius neighborhood analysis
----------------------------------

The ``nhood_analysis`` option constructs spatial neighborhoods directly from
cell coordinates.

Example:

.. code-block:: yaml

   nhood_analysis:
     nhood_radius_size: 50
     pixel_conversion_rate: 1
     max_distance_coocurrence: 200
     dpi: 150

The radius supplied to Squidpy is:

.. code-block:: text

   neighborhood_radius_size / pixel_conversion_rate

The analysis does not use PCA or expression-space nearest neighbors.

Neighborhood enrichment
~~~~~~~~~~~~~~~~~~~~~~~

Neighborhood enrichment compares the observed adjacency of major cell
populations with a permutation-based expectation.

Positive z-scores indicate population pairs found near each other more often
than expected; negative values indicate depletion/avoidance.

SpPrAn exports both the complete heatmap and population-highlighted views.

Interaction matrix
~~~~~~~~~~~~~~~~~~

The interaction matrix summarizes direct neighborhood relationships among major
cell populations using the fixed-radius spatial graph.

SpPrAn exports a complete matrix plus highlighted views for individual
populations.

Co-occurrence analysis
~~~~~~~~~~~~~~~~~~~~~~

Co-occurrence analysis describes how the spatial association of populations
changes with distance.

``max_distance_coocurrence`` controls the largest distance displayed on the
co-occurrence figure. It changes the displayed x-axis range rather than the
radius used to construct the neighborhood graph.

.. important::

   Neighborhood association is not equivalent to causal cell-cell
   communication. Spatial proximity may support a biological hypothesis, but
   mechanism requires additional evidence.

Advanced phenotype classification
---------------------------------

Set:

.. code-block:: yaml

   advanced_cell_classification: true

to further characterize discovered phenotypes using protein intensity.

The initial unsupervised phenotype remains defined by protein marker **positivity**.
Advanced classification then describes the intensity of positive protein markers.

For each sample and protein marker:

* protein-marker-negative cells are labeled ``Absent``;
* positive cells below the 25th percentile are ``Low``;
* positive cells from the 25th through 75th percentile are ``Medium``; and
* positive cells above the 75th percentile are ``High``.

Thresholds are calculated independently in each sample.

Outputs
~~~~~~~

Advanced classification generates:

``Protein marker intensity summary.csv``
   Low/medium/high counts and proportions for each positive protein marker within each
   discovered phenotype.

``Intensity profiles analysis.csv``
   Observed combinations of protein-marker-intensity states within each phenotype and
   their relative frequencies.

``Phenotype intensity summary.csv``
   Phenotype-level overview including the number of positive protein markers, number of
   observed intensity profiles, and frequency of the most common profile.

The intensity classes and threshold metadata are also stored in AnnData when
AnnData saving is enabled.

Biological interpretation
~~~~~~~~~~~~~~~~~~~~~~~~~

The low/medium/high categories are **sample-relative intensity descriptors**.
They should not be interpreted as universal expression thresholds across
patients or experiments without additional normalization/harmonization.

This distinction is important in tumor-microenvironment analyses: a cell can
belong to the same positivity-defined phenotype in two tumors while exhibiting
different relative protein-expression strengths.
