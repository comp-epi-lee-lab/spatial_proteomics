Installation
============

SpPrAn can be used as a packaged desktop application or installed from source.

Desktop application
-------------------

Pre-built Windows and macOS applications are distributed through GitHub
Releases. See :doc:`desktop_app`.

Install from source
-------------------

SpPrAn 0.3.1 requires Python 3.12 or newer.

.. code-block:: bash

   git clone https://github.com/comp-epi-lee-lab/spatial_proteomics.git
   cd spatial_proteomics
   python -m pip install -e .

Run:

.. code-block:: bash

   sppran --config config/config.yaml
