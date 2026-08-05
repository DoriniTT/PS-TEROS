============
Installation
============

Install psteros with the AiiDA QE dependencies:

.. code-block:: console

   $ python -m venv .venv
   $ . .venv/bin/activate
   $ python -m pip install -U pip
   $ python -m pip install -c constraints/aiida-qe-2026-08.txt \
       aiida-core aiida-workgraph aiida-quantumespresso aiida-pseudo pytest
   $ python -m pip install --no-deps -e .
   $ pytest -q tests/unit/test_public_api.py

``constraints/aiida-qe-2026-08.txt`` pins the tested AiiDA/WorkGraph socket
contract.  Use a dedicated environment: installing into an already-running
AiiDA profile can replace its workflow dependencies.

For calculations, configure an AiiDA profile, an aiida-quantumespresso ``pw.x``
code, and an aiida-pseudo family.  The maintained SnO2 campaign uses QE 7.6 on
the Bohr ``gpu_a100`` queue with ``SSSP/1.3/PBE/precision``.  VASP users may
install the optional ``.[vasp]`` extra and select ``backend="vasp"`` in the
same typed public workflow configuration; QE remains the primary path.
