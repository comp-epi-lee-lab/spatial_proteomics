Desktop application
===================

SpPrAn includes a Tkinter-based desktop application intended for users who
prefer not to work from the command line.

Distribution
------------

Packaged Windows and macOS applications are available through `GitHub
Releases <https://github.com/comp-epi-lee-lab/spatial_proteomics/releases>`_. Current supported platforms:

* **v0.2.5-beta.1**:

    * Windows 11 (x86_64 - 64-bit Intel or AMD processors)
    * Windows 11 (ARM64 / Copilot+ PCs)
    * macOS 26   (Apple Silicon: M1/M2/M3/M4)

Download and startup
--------------------

Download the ZIP file corresponding to your operating system and processor architecture from the **Assets** 
section in GitHub Releases.

Windows 11
^^^^^^^^^^

#. Download the appropriate Windows ZIP file:
    
    * Windows-x86_64 for 64-bit Intel or AMD processors
    * Windows-ARM64 for ARM-based Windows devices

#. Extract the contents of the ZIP file to a folder.
#. Open the extracted folder.
#. Double-click the SpPrAn GUI executable (.exe) to start the application.

Windows may display a security or SmartScreen warning because this beta build is not distributed through the 
Microsoft Store and may not be digitally signed. If this occurs, verify that you downloaded the application 
from this official GitHub Releases page before choosing to continue.

macOS
^^^^^

#. Download the macOS-arm64 ZIP file.
#. Extract the ZIP file.
#. Open the extracted folder.
#. Double-click the SpPrAn GUI application (.app) to start it.

Because this beta build may not yet be notarized or digitally signed by Apple, macOS may prevent the application 
from opening the first time.

If macOS reports that the application cannot be opened because the developer cannot be verified, you can open 
System Settings → Privacy & Security, locate the message referring to SpPrAn GUI, and choose Open Anyway.

Only bypass this warning if you downloaded the application directly from this project's official GitHub Releases 
page.

.. note::

    The first time it is started, the application takes about 1-2 minutes to fully load on any OS due to the 
    creation of cache data from the Matplotlib library. After that, the startup should take seconds.

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
