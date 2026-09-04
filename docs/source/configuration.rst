Configuration file
==================

SpPrAn is controlled by a YAML configuration file. Users can create it manually
from ``config/config_example.yaml`` or with the :doc:`config_builder`.

The configuration controls workspace paths, file format, protein markers,
primary cell types, hierarchical subtypes, colors, AnnData saving, overwrite
behavior, spatial plotting, custom plotting, and hierarchy orientation.

``protein_markers`` option is a list that contains the primary protein markers to be used in the 
cell type definitions. Protein markers name must have one of the following formats:
    * "Positivity - {``protein_marker_name``} (MV Cell)" 
    * "Positivity - {``protein_marker_name``}* (MV Cell)"
    * "Positivity - {``protein_marker_name``} (MV - NUC)" 
    * "Positivity - {``protein_marker_name``}* (MV - NUC)"
    * "Positivity - {``protein_marker_name``} (MV - CYTO)" 
    * "Positivity - {``protein_marker_name``}* (MV - CYTO)"

``cell_types`` option is a dictionary, where the keys are the cell type names and the values are 
a list with the detected intensity signal preserving the same order as ``protein_markers``. Primary 
phenotype definitions inside a cell type use ``1`` for positive, ``0`` for negative, and
``~``/``null`` for unconstrained markers.

``sub_cell_types`` option is a nested dictionary, where the first-level keys are the parent cell 
type and subtypes names and the values are dictionaries defining cell subtypes within that parent 
population. Each cell subtype name maps to a dictionary where the keys are protein markers and the 
values define the phenotype: ``1`` for positive, ``0`` for negative, and ``~``/``null`` for 
unconstrained markers. Cell type and subtype names cannot be the same.

``custom_colors`` option is a dictionary, where the keys are the cell type or subtype names 
and the values are hexadecimal color codes used to represent each population in plots or other 
visual outputs. Cell type and subtype names cannot be the same.

``spatial_plot`` option is a dictionary containing the parameters used to generate spatial plots. 
``dpi`` defines the plot resolution, ``scatter_point_size`` controls the size of the plotted cells 
or objects, and ``plot_options`` specifies the plotting mode.

The values for ``plot_options`` are ``overview``, ``highlight``, ``all``, ``none``, and ``custom``.
For more information, see :doc:`config_builder`. When ``plot_options`` is set to ``custom``, 
the ``custom_plot`` option is a list with the input files to plot, where each entry is a dictionary
that defines the input ``fileName`` and the corresponding plot ``type`` similar to ``plot_options``
(except for ``custom`` value).

An additional ``arrows_orientation`` option is available. For the cell hierarchy plot, 
``arrows_orientation`` will change the orientation of the arrows. ``contains`` value
will show arrows from up to down, and ``isContained`` value will show arrows from
down to up.
