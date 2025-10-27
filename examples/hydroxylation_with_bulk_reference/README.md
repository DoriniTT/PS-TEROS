# Hydroxylation Workflow with Bulk and Pristine Reference Calculations

This example demonstrates the **v2 hydroxylation workflow** with automatic bulk and pristine slab reference calculations for surface energy analysis.

## What's New in v2

- **Automatic bulk relaxation** from CIF file (ISIF=3 for cell optimization)
- **Automatic pristine slab relaxation** for reference surface energy (γ₀)
- **Four new outputs**: `bulk_structure`, `bulk_energy`, `pristine_structure`, `pristine_energy`
- Complete reference data for surface thermodynamics (Section S2)

## Requirements

- Relaxed slab structure (StructureData or VASP file)
- Bulk crystal structure (CIF file, StructureData, or PK)
- VASP code configured in AiiDA
- AiiDA daemon running

## Files

- `run_hydroxylation_v2.py`: Main workflow script with bulk reference
- `ag3po4.cif`: Bulk Ag₃PO₄ structure for reference calculations
- `README.md`: This file

## Usage

### 1. Prepare Structures

Ensure you have:
- Your relaxed slab structure (modify PK in script)
- Bulk structure CIF file (provided: `ag3po4.cif`)

### 2. Configure Parameters

Edit `run_hydroxylation_v2.py`:

```python
# Slab structure (choose one)
structure_pk = 1234  # Your slab PK
# OR
structure_file = Path("your_slab.vasp")

# Bulk structure (choose one)
bulk_cif_path='ag3po4.cif'
# OR
bulk_structure_pk=5678
# OR
bulk_structure=orm.load_node(5678)
```

### 3. Run Workflow

```bash
source ~/envs/aiida/bin/activate
cd /home/thiagotd/git/PS-TEROS/examples/hydroxylation_with_bulk_reference
python run_hydroxylation_v2.py
```

### 4. Monitor Progress

```bash
# Get workflow PK from output or workflow_pk.txt
verdi process show <PK>
verdi process report <PK>

# Watch live updates
watch -n 30 verdi process show <PK>
```

## Expected Workflow Steps

1. **Bulk relaxation** (ISIF=3): Cell + ionic optimization of ag3po4.cif
2. **Pristine slab relaxation** (ISIF=2): Ionic relaxation of unmodified slab
3. **Generate variants**: Create hydroxylated/vacancy structures
4. **Relax variants**: Parallel VASP relaxations

## Outputs

After completion, access results:

```python
from aiida import orm
from teros.core.surface_hydroxylation import organize_hydroxylation_results

node = orm.load_node(<workflow_pk>)
results = organize_hydroxylation_results(node)

# Reference calculations
ref = results['reference_data']
print(f"Bulk energy: {ref['bulk_energy']} eV")
print(f"Pristine energy: {ref['pristine_energy']} eV")
print(f"Bulk PK: {ref['bulk_structure_pk']}")
print(f"Pristine PK: {ref['pristine_structure_pk']}")

# Variant results
for r in results['successful_relaxations']:
    print(f"{r['name']}: {r['energy']:.6f} eV (coverage={r['coverage']:.2f})")
```

## Surface Energy Calculations

The reference data enables surface free energy calculations per Section S2:

- **E_bulk**: Use `ref['bulk_energy']`
- **γ₀**: Calculate from `ref['pristine_energy']` using equations 4-10
- **γ_modified**: Calculate for each variant using bulk and pristine references

## Comparison to v1

**v1 (run_hydroxylation_v1.py)**:
- Only relaxes modified surface structures
- No bulk or pristine reference calculations
- Cannot calculate surface energies

**v2 (run_hydroxylation_v2.py)**:
- Automatic bulk relaxation (ISIF=3)
- Automatic pristine slab relaxation
- Complete reference data for thermodynamics
- Ready for surface phase diagram generation

## Troubleshooting

**Error: "Bulk structure is required"**
- Provide one of: `bulk_cif_path`, `bulk_structure_pk`, or `bulk_structure`

**Error: "Could not read CIF file"**
- Check file exists: `ls ag3po4.cif`
- Check file format is valid CIF

**Bulk relaxation fails**
- Check `bulk_builder_inputs` has ISIF=3
- Verify VASP parameters are suitable for your system
- Check scheduler resources are adequate

## Related Documentation

- Design document: `/home/thiagotd/git/fosfato/calculos/hydroxylation/docs/plans/2025-10-27-bulk-pristine-reference-calculations-design.md`
- Section S2: `surface_energy_calc_procedure.tex`
- Main module: `teros/core/surface_hydroxylation/workgraph.py`
