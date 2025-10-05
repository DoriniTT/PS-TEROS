# Zone Approach Map-Zone Demonstration

This directory provides a concrete example of using the experimental
``Map`` context manager from **aiida-workgraph 1.0** to handle a
dynamic number of tasks. The workflow mirrors the PS-TEROS slab
relaxation use case: generate surface terminations on the fly and run a
VASP relaxation WorkChain for each one.

## Components

- ``workgraph.py`` – builds the production workgraph. It also contains a
  temporary shim that patches ``TaskManager.execute_map_task`` so that
  ``Map`` can iterate over dynamic task outputs until the upstream bug
  in aiida-workgraph is fixed.
- ``slabs_relax.py`` – CLI driver that loads the AiiDA profile and runs
  either the full VASP flow or a lightweight mock flow (``--mock``)
  that exercises the Map zone without launching real jobs.
- ``__init__.py`` – re-exports ``build_zone_workgraph`` and the mock
  builder for use elsewhere in PS-TEROS.

## Map Zone Flow

1. ``generate_slab_structures`` produces a dynamic namespace
   ``{'slabs': {label: StructureData}}`` that contains an unknown number
   of slabs.
2. On import, ``workgraph.py`` patches ``TaskManager.execute_map_task``
   so a dynamic namespace is correctly passed as the ``source`` input of
   the Map zone.
3. ``build_zone_workgraph`` loads the requested VASP code, creates a
   WorkChain task with ``task(WorkflowFactory('vasp.v2.vasp'))`` and maps
   it over the slab namespace. Each mapped iteration gathers the relaxed
   structure and the total energy (via ``extract_total_energy``) into the
   Map outputs.
4. The workgraph exposes ``generated_slabs``, ``relaxed_slabs`` and
   ``slab_energies`` as dynamic namespaces that can be inspected through
   AiiDA or consumed by downstream logic.

Because the relaxation step now submits the VASP WorkChain directly,
the provenance beneath the top-level ``WorkGraph`` consists of the
expected ``WorkChainNode<vasp.v2.vasp>`` children—no nested sub-
workgraphs.

## Running the Slab Workflow

```bash
source ~/envs/aiida/bin/activate
verdi profile set-default psteros
verdi daemon restart

cd teros/test_modules/zone_approach
python slabs_relax.py
```

The script loads ``structures/ag3po4.cif`` from the repository root and
launches one VASP WorkChain per generated slab via the Map zone. Inspect
the provenance with ``verdi process show <PK>``.

## Mock Mode (No VASP Required)

```bash
python slabs_relax.py --mock --mock-count 4 --mock-delta 0.5
```

The mock workflow generates simple ``Float`` values, maps over them, and
prints the gathered results. This is useful for testing the Map wiring
without relying on VASP inputs or cluster resources.

## Notes

- Remove the monkey patch in ``workgraph.py`` once aiida-workgraph ships
  a release where ``Map`` correctly propagates the ``source`` keyword for
  dynamic namespaces.
- The AiiDA daemon may warn that the RabbitMQ version is unsupported; see
  the linked wiki for upgrade instructions if you experience instability.
- Even the mock workflow creates transient AiiDA nodes. Ensure your user
  has write access to ``~/.aiida/access`` when running it.

