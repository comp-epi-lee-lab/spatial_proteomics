Release notes
=============

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
