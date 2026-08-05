.. _api:

=============
API Reference
=============

Public API
----------

``psteros`` deliberately exposes a small, typed public API.  QE through
``aiida-quantumespresso`` is the primary implementation; VASP configuration
is retained for established VASP studies.

.. automodule:: psteros
   :members:
   :undoc-members:

Workflow API
------------

.. automodule:: psteros.workflow
   :members:
   :undoc-members:

Configuration API
-----------------

.. automodule:: psteros.config
   :members:
   :undoc-members:

Compatibility API
-----------------

The pre-1.0 VASP builders are isolated under ``psteros.compat`` only for
existing projects.  New workflows must use the typed public API above;
retired interfaces have no supported import path.
