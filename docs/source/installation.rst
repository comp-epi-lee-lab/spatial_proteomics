Installation
============

SpPrAn can be used as a packaged desktop application or installed from source.

Desktop application
-------------------

Precompiled Windows and macOS builds are distributed through GitHub Releases
when available. See :doc:`desktop_app`.

Install from source
-------------------

SpPrAn requires Python 3.12 or newer.

.. code-block:: bash

   git clone https://github.com/comp-epi-lee-lab/spatial_proteomics.git
   cd spatial_proteomics
   python -m pip install -e .

Run:

.. code-block:: bash

   sppran --config config/config.yaml
