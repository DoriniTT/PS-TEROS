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

* A supported DFT code (VASP, Quantum ESPRESSO, or CP2K)
* Appropriate AiiDA plugins for these codes

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

The sources for AiiDA-TEROS can be downloaded from the `Github repo`_.

.. code-block:: console

    $ git clone git@github.com:DoriniTT/aiida-teros.git
    $ cd aiida-teros
    $ pip install -e .

This will install the package in development mode, allowing you to modify the code and have changes immediately available.

.. _Github repo: https://github.com/DoriniTT/aiida-teros
