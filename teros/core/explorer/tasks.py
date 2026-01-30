"""Calcfunction tasks for the explorer module."""

import re
import os
import glob
import shutil
import tempfile
import subprocess
from pathlib import Path

from aiida import orm
from aiida_workgraph import task


@task.calcfunction
def compute_dynamics(
    structure: orm.StructureData,
    fix_type: orm.Str,
    fix_thickness: orm.Float,
    fix_elements: orm.List = None,
) -> orm.Dict:
    """
    Compute selective dynamics (positions_dof) for a structure.

    This calcfunction is used when the structure is not known at build time
    (e.g., from a previous stage's output), so we need to compute the fixed
    atoms at runtime.

    Args:
        structure: Input structure to analyze
        fix_type: Where to fix atoms ('bottom', 'center', 'top')
        fix_thickness: Thickness in Angstroms for fixing region
        fix_elements: Optional list of element symbols to restrict fixing to

    Returns:
        Dict with 'positions_dof' array for VaspWorkChain dynamics input
    """
    from teros.core.fixed_atoms import get_fixed_atoms_list

    fix_type_val = fix_type.value
    fix_thickness_val = fix_thickness.value
    fix_elements_val = fix_elements.get_list() if fix_elements is not None else None

    # Get list of atom indices to fix (1-based)
    fixed_atoms_list = get_fixed_atoms_list(
        structure=structure,
        fix_type=fix_type_val,
        fix_thickness=fix_thickness_val,
        fix_elements=fix_elements_val,
    )

    # Create positions_dof array: True = relax, False = fix
    num_atoms = len(structure.sites)
    positions_dof = []

    for i in range(1, num_atoms + 1):  # 1-based indexing
        if i in fixed_atoms_list:
            positions_dof.append([False, False, False])  # Fix atom
        else:
            positions_dof.append([True, True, True])  # Relax atom

    return orm.Dict(dict={'positions_dof': positions_dof})


@task.calcfunction
def extract_energy(misc: orm.Dict, retrieved: orm.FolderData = None) -> orm.Float:
    """
    Extract total energy from VASP misc output or OUTCAR.

    Args:
        misc: VASP misc output Dict containing energy data
        retrieved: VASP retrieved FolderData (optional) to parse OUTCAR if misc fails

    Returns:
        Total energy as Float (eV)
    """
    misc_dict = misc.get_dict()

    # Navigate to total_energies if present
    energy_dict = misc_dict
    if 'total_energies' in misc_dict:
        energy_dict = misc_dict['total_energies']

    # Try multiple keys in order of preference
    for key in ('energy_extrapolated', 'energy_no_entropy', 'energy'):
        if key in energy_dict:
            return orm.Float(float(energy_dict[key]))

    # If no recognized key found in misc, try to parse from retrieved OUTCAR
    if retrieved is not None:
        try:
            content = retrieved.get_object_content('OUTCAR')
            # Look for "free  energy   TOTEN  =       -832.63657516 eV"
            matches = re.findall(r'free\s+energy\s+TOTEN\s+=\s+([-\d.]+)', content)
            if matches:
                return orm.Float(float(matches[-1]))
        except Exception:
            pass

    # If no recognized key found, raise error with available keys
    available = ', '.join(sorted(energy_dict.keys()))
    raise ValueError(f'Unable to find total energy in misc output or OUTCAR. Available keys: {available}')


@task.calcfunction(outputs=['charges', 'acf', 'bcf', 'avf'])
def run_bader_analysis(retrieved: orm.FolderData, structure: orm.StructureData) -> dict:
    """
    Run Bader charge analysis on AECCAR files from a VASP SCF calculation.

    This calcfunction:
    1. Extracts AECCAR0, AECCAR2, CHGCAR from the retrieved FolderData
    2. Sums AECCAR0 + AECCAR2 using pymatgen to create CHGCAR_sum
    3. Runs the bader binary: ``bader CHGCAR -ref CHGCAR_sum``
    4. Parses ACF.dat and returns charges + all .dat files

    Args:
        retrieved: FolderData from a VASP SCF calculation that produced
                   AECCAR0, AECCAR2, and CHGCAR (requires ``laechg: True``
                   in INCAR and these files in ADDITIONAL_RETRIEVE_LIST)
        structure: StructureData with the same atom ordering as the VASP
                   calculation. Used to map atoms to elements and look up
                   valence electron counts (ZVAL) from the OUTCAR.

    Returns:
        dict with:
            - 'charges': orm.Dict with parsed ACF.dat data (atoms, charges, volumes,
              element, valence, bader_charge per atom)
            - 'acf': orm.SinglefileData for enriched ACF.dat (with ELEMENT, VALENCE,
              BADER_CHARGE columns)
            - 'bcf': orm.SinglefileData for BCF.dat
            - 'avf': orm.SinglefileData for AVF.dat
            - Additional .dat files as SinglefileData with lowercase keys
    """
    from pymatgen.io.vasp.outputs import Chgcar

    # Find the bader binary
    bader_path = Path.home() / '.local' / 'bin' / 'bader'
    if not bader_path.exists():
        raise FileNotFoundError(
            f"Bader binary not found at {bader_path}. "
            f"Install it to ~/.local/bin/bader"
        )

    # Create a temporary working directory
    tmpdir = tempfile.mkdtemp(prefix='bader_')
    try:
        # Extract required files from retrieved FolderData
        required_files = ['AECCAR0', 'AECCAR2', 'CHGCAR', 'OUTCAR']
        for fname in required_files:
            try:
                content = retrieved.get_object_content(fname, mode='rb')
            except (FileNotFoundError, OSError):
                raise FileNotFoundError(
                    f"File '{fname}' not found in retrieved FolderData (PK {retrieved.pk}). "
                    f"Ensure the SCF stage has 'laechg': True in INCAR and "
                    f"['AECCAR0', 'AECCAR2', 'CHGCAR'] in the retrieve list."
                )
            filepath = os.path.join(tmpdir, fname)
            with open(filepath, 'wb') as f:
                f.write(content)

        # Sum AECCAR0 + AECCAR2 → CHGCAR_sum using pymatgen
        aeccar0 = Chgcar.from_file(os.path.join(tmpdir, 'AECCAR0'))
        aeccar2 = Chgcar.from_file(os.path.join(tmpdir, 'AECCAR2'))
        chgcar_sum = aeccar0 + aeccar2
        chgcar_sum_path = os.path.join(tmpdir, 'CHGCAR_sum')
        chgcar_sum.write_file(chgcar_sum_path)

        # Run bader binary
        result = subprocess.run(
            [str(bader_path), 'CHGCAR', '-ref', 'CHGCAR_sum'],
            cwd=tmpdir,
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"Bader analysis failed (return code {result.returncode}).\n"
                f"stdout: {result.stdout}\n"
                f"stderr: {result.stderr}"
            )

        # Parse ZVAL from OUTCAR using pymatgen
        from pymatgen.io.vasp import Outcar
        zval_dict = {}
        try:
            outcar = Outcar(os.path.join(tmpdir, 'OUTCAR'))
            outcar.read_pseudo_zval()
            zval_dict = outcar.zval_dict  # e.g. {'Sn': 14.0, 'O': 6.0}
        except Exception:
            pass  # If OUTCAR parsing fails, proceed without ZVAL

        # Build element list from structure (matches ACF.dat atom order)
        elements = [site.kind_name for site in structure.sites]

        # Collect all .dat files
        dat_files = glob.glob(os.path.join(tmpdir, '*.dat'))

        # Parse ACF.dat for charges
        acf_path = os.path.join(tmpdir, 'ACF.dat')
        charges_data = _parse_acf_dat(acf_path)

        # Enrich charges data with element, valence, and bader_charge
        for idx, atom in enumerate(charges_data.get('atoms', [])):
            if idx < len(elements):
                elem = elements[idx]
                atom['element'] = elem
                valence = zval_dict.get(elem, 0.0)
                atom['valence'] = valence
                atom['bader_charge'] = valence - atom['charge']

        # Enrich ACF.dat file with ELEMENT, VALENCE, BADER_CHARGE columns
        _enrich_acf_dat(acf_path, elements, zval_dict)

        # Build output dict
        outputs = {}
        outputs['charges'] = orm.Dict(dict=charges_data)

        # Store each .dat file as SinglefileData
        for dat_path in dat_files:
            basename = os.path.basename(dat_path)
            # Key: lowercase name without extension (acf, bcf, avf, etc.)
            key = basename.replace('.dat', '').lower()
            outputs[key] = orm.SinglefileData(file=dat_path)

        return outputs

    finally:
        # Clean up temporary directory
        shutil.rmtree(tmpdir, ignore_errors=True)


def _parse_acf_dat(filepath: str) -> dict:
    """
    Parse ACF.dat file from Bader analysis.

    ACF.dat format:
        Line 1: header
        Line 2: separator (dashes)
        Lines 3..N+2: atom data (index, x, y, z, charge, min_dist, volume)
        Line N+3: separator (dashes)
        Lines after: summary

    Args:
        filepath: Path to ACF.dat file

    Returns:
        dict with:
            - atoms: list of dicts with keys: index, x, y, z, charge, min_dist, volume
            - total_charge: float (sum of all Bader charges)
            - vacuum_charge: float (charge not assigned to any atom)
            - vacuum_volume: float (volume not assigned to any atom)
    """
    atoms = []
    total_charge = 0.0
    vacuum_charge = 0.0
    vacuum_volume = 0.0

    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Skip header (line 0) and separator (line 1)
    for line in lines[2:]:
        line = line.strip()
        if not line or line.startswith('-'):
            continue

        # Check for summary lines
        parts = line.split()

        # Summary lines: "Number of electrons:" or "Vacuum charge:" etc.
        if 'electrons' in line.lower():
            # "Number of electrons:   192.00000"
            try:
                total_charge = float(parts[-1])
            except (ValueError, IndexError):
                pass
            continue
        if 'vacuum' in line.lower() and 'charge' in line.lower():
            try:
                vacuum_charge = float(parts[-1])
            except (ValueError, IndexError):
                pass
            continue
        if 'vacuum' in line.lower() and 'volume' in line.lower():
            try:
                vacuum_volume = float(parts[-1])
            except (ValueError, IndexError):
                pass
            continue

        # Atom data line: index x y z charge min_dist volume
        try:
            if len(parts) >= 7:
                atom = {
                    'index': int(parts[0]),
                    'x': float(parts[1]),
                    'y': float(parts[2]),
                    'z': float(parts[3]),
                    'charge': float(parts[4]),
                    'min_dist': float(parts[5]),
                    'volume': float(parts[6]),
                }
                atoms.append(atom)
        except (ValueError, IndexError):
            continue

    return {
        'atoms': atoms,
        'total_charge': total_charge,
        'vacuum_charge': vacuum_charge,
        'vacuum_volume': vacuum_volume,
    }


def _enrich_acf_dat(filepath: str, elements: list, zval_dict: dict) -> None:
    """
    Rewrite ACF.dat in place with added ELEMENT, VALENCE, and BADER_CHARGE columns.

    Args:
        filepath: Path to ACF.dat file
        elements: List of element symbols matching atom order
        zval_dict: Dict mapping element symbol to ZVAL (valence electrons)
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    if len(lines) < 3:
        return

    new_lines = []
    atom_idx = 0

    for i, line in enumerate(lines):
        stripped = line.strip()

        if i == 0:
            # Header line: append new column headers
            new_lines.append(
                stripped + '  ELEMENT  VALENCE  BADER_CHARGE\n'
            )
        elif stripped.startswith('-'):
            # Separator line: extend dashes
            new_lines.append(
                stripped + '----------------------------\n'
            )
        elif (
            'electrons' in stripped.lower()
            or ('vacuum' in stripped.lower() and 'charge' in stripped.lower())
            or ('vacuum' in stripped.lower() and 'volume' in stripped.lower())
            or not stripped
        ):
            # Summary/footer lines: keep as-is
            new_lines.append(line)
        else:
            # Atom data line
            parts = stripped.split()
            try:
                if len(parts) >= 7 and atom_idx < len(elements):
                    elem = elements[atom_idx]
                    valence = zval_dict.get(elem, 0.0)
                    charge = float(parts[4])
                    bader_charge = valence - charge
                    new_lines.append(
                        f'{stripped}  {elem:>7s}  {valence:7.1f}  '
                        f'{bader_charge:12.6f}\n'
                    )
                    atom_idx += 1
                else:
                    new_lines.append(line)
            except (ValueError, IndexError):
                new_lines.append(line)

    with open(filepath, 'w') as f:
        f.writelines(new_lines)
