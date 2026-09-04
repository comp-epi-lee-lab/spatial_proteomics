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

Typical workflow
----------------

#. Export object-level measurements from Visiopharm®.
#. Create or update a YAML configuration.
   #. Define biologically appropriate marker rules.
   #. Add hierarchical subtypes where supported by the marker panel.
   #. Select spatial plotting/output options.
#. Run SpPrAn from the desktop application or command line.
#. Review population counts, hierarchy reports, and spatial maps.
