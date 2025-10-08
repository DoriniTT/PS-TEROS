
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
