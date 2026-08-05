============
Installation
============

Start from a source checkout and a fresh virtual environment. The steps below
install psteros and its declared Python dependencies; they do not configure a
remote computer or submit a calculation.

Before you begin
----------------

You need Python 3.11 or newer and a checkout of the psteros repository. Use a
dedicated virtual environment so that project dependencies do not unexpectedly
change an AiiDA environment you already use for other work.

Create the environment
----------------------

From the repository root, create and activate a virtual environment:

.. code-block:: console

   $ python -m venv .venv
   $ . .venv/bin/activate
   $ python -m pip install --upgrade pip

Install psteros and its declared dependencies:

.. code-block:: console

   $ python -m pip install .

Verify the installation
-----------------------

Check that Python can import the package and report its version:

.. code-block:: console

   $ python -c "import psteros; print(psteros.__version__)"
   1.0.0

Prepare AiiDA before using a calculation graph
-----------------------------------------------

The tutorial builds a graph locally, but a real calculation also needs an AiiDA
profile with three project-specific pieces:

* a configured **computer**, which represents the machine or cluster that will
  run the calculation;
* a registered Quantum ESPRESSO ``quantumespresso.pw`` **code** on that
  computer; and
* an ``aiida-pseudo`` **pseudopotential family** compatible with the method you
  plan to use.

AiiDA stores those choices as part of the calculation record. Follow AiiDA's
`guide to configuring and running external codes`_ for the computer and code
setup. Use the identifiers from *your* profile in the psteros examples; do not
copy a label, queue, or resource request from another project as a default.

.. _guide to configuring and running external codes:
   https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/run_codes.html

Optional VASP support
---------------------

The main supported workflow uses Quantum ESPRESSO. If you are continuing an
existing VASP study, install its optional adapter in the same environment:

.. code-block:: console

   $ python -m pip install '.[vasp]'

What next?
----------

Once the package imports and your AiiDA profile is ready, continue with the
:doc:`first tutorial <tutorial>`. It builds one graph without submitting work,
so you can inspect the recipe before requesting compute time.
