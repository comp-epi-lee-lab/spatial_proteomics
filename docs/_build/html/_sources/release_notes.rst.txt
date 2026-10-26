Release notes
=============

SpPrAn 0.3.1
------------

Version 0.3.1 expands the unsupervised/full workflow introduced in v0.3.0.

Added
~~~~~

* Continuous protein-marker intensity spatial plots.
* Fixed-radius spatial-neighborhood construction.
* Neighborhood-enrichment overview and highlighted plots.
* Cell-population interaction matrices and highlighted views.
* Distance-dependent co-occurrence plots.
* Advanced phenotype classification using sample-relative low/medium/high
  protein-intensity categories.
* Aggregate intensity-profile summaries for discovered phenotypes.
* Compact phenotype IDs and mapping tables for phenotype-level intensity output
  filenames.

Configuration Builder
~~~~~~~~~~~~~~~~~~~~~

A new **Optional spatial and phenotype analyses** section exposes controls for:

* intensity-spatial plot DPI and scatter-point size;
* neighborhood radius;
* pixel/coordinate conversion rate;
* maximum co-occurrence display distance;
* neighborhood plot DPI; and
* advanced phenotype classification.

Interpretation
~~~~~~~~~~~~~~

The new features intentionally separate three analytical levels:

* marker positivity defines phenotype identity;
* continuous/categorical intensity describes expression strength; and
* neighborhood analyses describe spatial organization.

Backward compatibility
~~~~~~~~~~~~~~~~~~~~~~

The new v0.3.1 settings are optional. Existing v0.3.0 configurations that omit
them do not enable the additional analyses.

SpPrAn 0.3.0
------------

Version 0.3.0 introduces a mode-based analysis architecture.

Major additions include:

* supervised, unsupervised, and full analysis dispatch;
* exploratory marker-pattern discovery within selected major cell types;
* cell-type proportion matrices and optional pie plots;
* phenotype-proportion matrices;
* optional two-group phenotype comparison;
* Mann--Whitney U and Welch t-test summaries with FDR correction;
* optional phenotype-specific protein-marker intensity summaries;
* full phenotype discovery that converts discovered marker patterns into
  subtype definitions and continues through the supervised workflow; and
* expanded Configuration Builder controls.

Backward compatibility
~~~~~~~~~~~~~~~~~~~~~~

Configurations that predate ``type_of_analysis`` are interpreted as supervised
configurations.
