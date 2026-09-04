Desktop application
===================

SpPrAn includes a Tkinter-based desktop application intended for users who
prefer not to work from the command line.

Distribution
------------

Packaged Windows and macOS applications are available through GitHub
Releases. Current suported platforms:

* Windows 11 (x86_64)
* Windows 11 (ARM64)
* macOS 26 (ARM64)

Main interface
--------------

The desktop application provides controls to:

* create, edit, or open a configuration file;
* choose input and output directories;
* inspect the selected configuration;
* run the SpPrAn analysis; and
* view progress and error messages.

Configuration Builder
---------------------

Selecting **Create a new or edit an existing configuration file** opens the
bundled :doc:`config_builder` in the user's default browser.

Running an analysis
-------------------

Before an analysis can be launched, the desktop application requires a
configuration file, an input directory, and an output directory. It then calls
the same underlying SpPrAn analysis pipeline used by the command-line
interface.
