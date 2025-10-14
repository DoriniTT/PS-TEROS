.. _installation:

============
Installation
============

Requirements
------------

TEROS requires:

* Python 3.9 or higher
* AiiDA Core
* AiiDA WorkGraph
* pymatgen
* numpy

For running actual calculations, you will also need:

* VASP (currently the only supported DFT code)
* AiiDA-VASP plugin

.. note::
   **Current DFT Code Support**: The current version only supports VASP. CP2K and Quantum ESPRESSO are planned as future updates. If you need PS-TEROS code for these versions, please use the legacy version in the ``legacy/`` directory.

AiiDA Installation
------------------

The full and complete AiiDA installation process can be found here:
https://aiida.readthedocs.io/projects/aiida-core/en/stable/installation/guide_complete.html

For an easy setup, you can follow the quick installation guide:
https://aiida.readthedocs.io/projects/aiida-core/en/stable/installation/guide_quick.html

This typically involves:

.. code-block:: console

    $ pip install aiida-core
    $ verdi presto
    $ verdi status

From Source
-----------

The sources for PS-TEROS can be downloaded from the `Github repo`_.

.. code-block:: console

    $ git clone git@github.com:DoriniTT/PS-TEROS.git
    $ cd PS-TEROS
    $ pip install -e .

This will install the package in development mode, allowing you to modify the code and have changes immediately available.

.. _Github repo: https://github.com/DoriniTT/PS-TEROS
