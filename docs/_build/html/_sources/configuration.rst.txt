Configuration file
==================

SpPrAn is controlled by a YAML configuration file. Users can create it manually from ``config/config_example.yaml``
or with the :doc:`config_builder`.

Analysis type
-------------

``type_of_analysis`` accepts:

``supervised``
   User-defined major cell types and hierarchical subtypes.

``unsupervised``
   User-defined major cell types followed by exploratory marker-pattern
   discovery within selected parents.

``full``
   Unsupervised discovery followed by conversion of discovered patterns into
   subtype definitions and execution of the supervised workflow.

Missing or invalid values default to ``supervised`` for backward compatibility.

Shared settings
---------------

``workspace``
   Input directory, output directory, and ``tsv``/``csv`` file type.

``protein_markers``
   List that contains the primary protein markers to be used in the cell type definitions. Protein 
   markers name must have one of the following formats:
        * "Positivity - {``protein_marker_name``} (MV Cell)" 
        * "Positivity - {``protein_marker_name``}* (MV Cell)"
        * "Positivity - {``protein_marker_name``} (MV - NUC)" 
        * "Positivity - {``protein_marker_name``}* (MV - NUC)"
        * "Positivity - {``protein_marker_name``} (MV - CYTO)" 
        * "Positivity - {``protein_marker_name``}* (MV - CYTO)"

``cell_types``
   Dictionary, where the keys are the cell type names and the values are 
   a list with the detected intensity signal preserving the same order as ``protein_markers``. Primary 
   phenotype definitions inside a cell type use ``1`` for positive, ``0`` for negative, and
   ``~``/``null`` for unconstrained markers.

``custom_colors``
   Dictionary where the keys are the cell type or subtype names 
   and the values are hexadecimal color codes used to represent each population in plots or other 
   visual outputs. Cell type and subtype names cannot be the same.

``locally_save_anndata_files``
   Whether annotated AnnData objects should be saved. ``true`` or ``false`` values are only accepted.

``overwrite_existing_files``
   Whether existing output files may be replaced.``true`` or ``false`` values are only accepted.

``arrows_orientation``
   For the cell hierarchy plot, ``arrows_orientation`` will change the orientation of the arrows. 
   ``contains`` value will show arrows from up to down, and ``isContained`` value will show arrows 
   from down to up.

Supervised-only settings
------------------------

``sub_cell_types``
   nested dictionary, where the first-level keys are the parent cell 
   type and subtypes names and the values are dictionaries defining cell subtypes within that parent 
   population. Each cell subtype name maps to a dictionary where the keys are protein markers and the 
   values define the phenotype: ``1`` for positive, ``0`` for negative, and ``~``/``null`` for 
   unconstrained markers. Cell type and subtype names cannot be the same.

``spatial_plot``
   dictionary containing the parameters used to generate spatial plots. 
   ``dpi`` defines the plot resolution, ``scatter_point_size`` controls the size of the plotted cells 
   or objects, and ``plot_options`` specifies the plotting mode.
   
   The values for ``plot_options`` are ``overview``, ``highlight``, ``all``, ``none``, and ``custom``.
   For more information, see :doc:`config_builder`. When ``plot_options`` is set to ``custom``, 
   the ``custom_plot`` option is a list with the input files to plot, where each entry is a dictionary
   that defines the input ``fileName`` and the corresponding plot ``type`` similar to ``plot_options``
   (except for ``custom`` value).

Unsupervised/full settings
--------------------------

``cell_types_to_explore``
   List with the cell types within which phenotypes should be discovered. Cell types names 
   must be exactly the same as in ``cell_types``, but not necessarily in the same order.

``cell_types_to_explore_short_name``
   List with the short names corresponding to the selected cell types in ``cell_types_to_explore``. Short names 
   must be different each other.

``exploratory_protein_markers``
   List with the protein markers used for exploratory phenotype discovery. Protein marker names 
   must be exactly the same as appear in the input files.

``exploratory_protein_markers_short_name``
   List with the short names corresponding to the exploratory protein markers. Short names 
   must be different each other.

``want_to_plot_pie_plots``
   Boolean controlling sample-level composition pie plots. ``true`` or ``false`` values are only accepted.

``dict_ctrl`` and ``dict_trtmnt``
   Optional filename-to-group dictionary mappings for two-group comparison. Keys are filenames in the input folder 
   and values is the corresponding group name. Group names must be different between them, and each dictionary
   must not contain a similar key. Both dictionaries can be empty only at the same time.

``protein_markers_for_intensity``
   Optional list with the protein markers for summarized intensity measurements. Protein marker names 
   must be exactly the same as appear in the input files. This list can be empty.

Validation
----------

SpPrAn validates mode-specific settings before analysis.

See :doc:`analysis_modes` for analytical interpretation.
