SpPrAn Configuration Builder
============================

The **SpPrAn Configuration Builder** is the local browser-based interface in
``config_builder.html``. It can be opened directly or launched from the desktop
application.

.. important::

   Files selected in the Configuration Builder are not uploaded to an external
   service. Representative data files are read locally to identify column
   names, and folders selected are used locally to obtain filenames.

Choose an analysis mode
-----------------------

SpPrAn begins with a mode-selection view.

**Supervised mode**
   Users define cell types and subtypes directly from protein-marker rules.

**Unsupervised mode**
   Users define major cell types using protein-marker rules, while subtypes are
   discovered from observed protein-marker patterns.

See :doc:`analysis_modes` for the analytical distinction.

Shared setup
------------

Both workflows require:

#. a representative ``.tsv`` or ``.csv`` file, or an existing YAML
   configuration;
#. workspace input/output settings; and
#. primary protein markers.


Settings 1. Restore or start
----------------------------

Loading an existing ``.yaml``/``.yml`` configuration allows users to edit an existing one or take it as a 
template for a new one.

Selecting a representative ``.tsv``/``.csv`` file allows users to start a new configuration file identifying 
available positivity columns in their input data.

.. important::

   * Configuration files from before 0.3.0 are treated as **supervised mode configuration files** for backward compatibility when restored.
   * Configuration files from before 0.3.0 can be used in the unsupervised mode restoring option. However, it will restore partial information.

Settings 2. Workspace
---------------------

Specify the input directory, output directory, and common input format. For desktop
application users, input and output directory can be selected through the corresponding
app.

.. important::

   * All the files in the input directory must be the same type. If not, SpPrAn will ignore those different from the selected input-file format.
   * Input file names must end with the suffix "**_objects.{``filetype``}**" to be analyzed by SpPrAn. Otherwise, code will not continue with the analysis.

.. note::

   We strongly recommend to only have the files to be analized inside the input directory and no other
   non-intented ``.tsv``/``.csv`` file, even if they are inside a folder.

.. note::

   We strongly recommend naming the input files with the suffix "**XY_objects.{``filetype``}**", 
   where **X** and **Y** are together the "Sample Id" in the analysis. **X** and **Y** can be any valid 
   character allowed by the OS. Because "**XY**" is the Id of the sample, each file must contain 
   different Ids.

Settings 3. Protein markers
---------------------------

Select positivity columns obtained via the representative file to define major 
populations (cell types). Other protein markers will remain available for other definitions 
if needed.

.. important::

   The allowed formats for columns in input data to be accepted as "positive columns" are:

   * "Positivity - {``protein_marker_name``} (MV Cell)" 
   * "Positivity - {``protein_marker_name``}* (MV Cell)"
   * "Positivity - {``protein_marker_name``} (MV - NUC)" 
   * "Positivity - {``protein_marker_name``}* (MV - NUC)"
   * "Positivity - {``protein_marker_name``} (MV - CYTO)" 
   * "Positivity - {``protein_marker_name``}* (MV - CYTO)"

.. important::

   When loading an existing configuration file, only protein marker used in the
   file will appear. Please select a representative file to retrieve and make 
   available all the protein markers you may need.

.. note::

   **In unsupervised configuration, user must define a short name for protein markers.** The Configuration Builder
   offers a default name that user can modify. Protein marker short names cannot have the same name between them.

Supervised configuration
------------------------

In supervised mode, users define:

* primary cell types and optional hierarchical subtypes using additional protein markers;
* AnnData saving, overwrite behavior and spatial plotting options.


"Cell types and subtypes" section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this settings section, users define each major cell type using protein markers 
that are **Present**, **Absent**, or **Not applicable**, then optionally create 
nested subtypes using additional markers.

Protein markers selected in settings option 3 will be used to define cell types. Please 
select introduce a name for the cell type and select the detected intensity signal for the protein 
markers to define it.

.. note::
   
   * Cell types and subtypes' name can be modified clicking in **Edit name** button.
   * Cell types and subtypes cannot have the same name or same detected intensity signal for the protein marker definition.

Cell subtypes can be added to each defined cell type by checking the **Include cell subtypes** 
checkbox. Cell subtypes will share the same protein markers used for defining the parent cell 
type. After naming a cell subtype, remaining protein markers can be selected for defining this subtype.

To define a cell subtype from a cell subtype, follow the steps above.

.. note::

   A parent-child relationship represents a nested marker-based phenotype
   definition, not automatically a developmental lineage. For example:

   .. code-block:: bash
   
      Immune cells
      └── T cells
         └── CD8 T cells

   represents nested marker-based phenotype definitions, not by itself a developmental lineage 
   or differentiation trajectory.

Click the **Collapse** button to temporarily hide the visualization of the cell type or 
subtype in the webpage. This is especially useful when there are many defined cell 
types and subtypes.

*Other cells* will collect all the cells that does not belong to any cell type or subtype. In
each parent-child relationship, it can appear *Other* {``parent_name``} corresponding to the cells
that belong to the parent but no to any of the children.

"Output and spatial plotting" section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this settings section, users may configure output overwriting, AnnData saving, plot type, DPI, 
point size, and optional per-sample custom plotting.

* For **Overwrite existing outputs?** option, if *Yes* is selected, previous generated output files in the output path will be replaced by the new ones. If *No* is selected, SpPrAn will skip generating all the output files that already exist, just creating the non-existent.
* The AnnData files can be used for generating plots inside and outside this tool, and follow-up analyses outside this tool using Python or R. For **Save AnnData files?** option, if *Yes* is selected, AnnData files will be saved in the output path.
* There are 5 options for plotting spatial plots:
   * **Overview**: SpPrAn will generate the spatial plots considering all annotations together: in 'General' all the cell types in the spatial plot, and all the children for each parent.
   * **Highlight**: SpPrAn will generate the spatial plots considering a single annotation vs all other cells: if "Immune cells" is a cell type, it will be Immune cells against all non-immune cells.
   * **Both**: SpPrAn will generate the spatial plots considering the overview and highlight options at the same time. This is the default option.
   * **No spatial plotting**: SpPrAn won't generate spatial plots for all samples in the input directory.
   * **Custom plotting**: SpPrAn will generate custom spatial plot types for each sample in the input directory.

For custom plotting, the builder reads matching filenames from a selected local
folder (input folder) and lets the user assign each file an overview, highlight, both, or no
plot. After assigning plot types for each file, click **Save custom plotting options** 
button for effectively save the assignations.

.. important::

   This is the only option that needs to be saved through clicking a button; every other 
   setting update is automatically saved.

.. important::

   Verify that the corrected input file format is selected before selecting the input folder. If not, you might see no files.

Unsupervised configuration
--------------------------

After defining major cell types, the unsupervised view asks which populations
should be explored and which protein markers should define exploratory phenotypes. Discovered phenotypes represent 
observed combinations of these protein marker states within the selected parent populations.

"Cell types to explore" section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this settings section, users define each major cell type using protein markers 
that are **Present**, **Absent**, or **Not applicable**. Cell types also need a short name, that it is used 
in result and annotation labels.

.. note::

   It is not possible to define nested cell subtypes in this mode.

"Phenotype discovery" section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this settings section, users choose how broadly SpPrAn should search for protein-marker-defined phenotypes
 within the defined cell types.

The **Full phenotype discovery** option, if checked, discovers protein marker patterns, stores them as subtype definitions, and
continues with the supervised workflow. If not checked, SpPrAn will do only the protein-marker-defined phenotype discovery.

The **Plot cell-type proportion pie plots?** option controls whether sample-level composition pie plots
are generated together with the cell-type proportion matrices.

The **Overwrite existing outputs?** and **Save AnnData files?** are exactly the same as in the supervised configuration.

"Optional spatial and phenotype analyses" section (**New in v0.3.1**)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this settings section, users optionally may choose additional analyses and plots to generate after phenotype discovery. 

The **Protein marker intensity spatial plots** option, if checked, creates one continuous-intensity spatial map per sample and
protein marker.

The **Neighborhood analysis** option, if checked, measures which cell populations tend to be near each other, how frequently they 
interact, and how their spatial association changes across distance.

.. note::

   The pixel conversion rate can be interpreted as 1:"Pixel conversion rate". That is because the neighborhood radius supplied 
   to the code is "Neighborhood radius size" / "Pixel conversion rate". Therefore, if the input data unit is 1:10, consider 
   modify "Pixel conversion rate" to 10 in order to correctly measure the neighborhood radius.

   The pixel conversion rate also have an impact in the maximum distance for co-occurrence plot.

The **Advanced phenotype classification** option, if checked, characterizes previously discovered phenotypes using low, medium, and 
high protein-intensity levels within each sample. This adds expression-strength information to discovered phenotypes without
changing the original positivity-based phenotype definition.

"Group-comparison analysis" section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this settings section, users optionally assign samples to two groups to compare discovered phenotypes between 
them. Here, the builder reads matching filenames from a selected local folder (input folder) and lets the user 
assign each file to each group. 

.. important::

   * If comparison is enabled, both groups must contain samples and use different group labels.
   * Samples may also remain unassigned.
   * Verify that the corrected input file format is selected before selecting the input folder.

"Intensity box-plot data" section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this settings section, users may select markers whose intensity measurements should be summarized
within discovered phenotypes for external downstream boxplot analysis.

"Configuration preview and download" section
--------------------------------------------

Preview the generated configuration in the browser. Download it for use with either
the desktop application or in the command line:

.. code-block:: bash

   sppran --config path/to/config.yaml
