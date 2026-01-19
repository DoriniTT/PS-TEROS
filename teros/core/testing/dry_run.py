"""Dry-run utilities for testing without cluster submission.

This module provides utilities to generate VASP input files without
actually submitting jobs to a cluster. Uses AiiDA's built-in dry_run mode.

Example:
    >>> from teros.core.testing import dry_run_vasp
    >>> from pymatgen.core import Structure, Lattice
    >>> structure = Structure(Lattice.cubic(4.0), ['Si', 'Si'],
    ...                       [[0, 0, 0], [0.25, 0.25, 0.25]])
    >>> result = dry_run_vasp(structure, builder_inputs, 'vasp@localhost')
    >>> result.print_incar()  # Inspect generated INCAR
"""

from pathlib import Path
from dataclasses import dataclass
from typing import Optional, Union, Dict, Any
import tempfile
import shutil
import os


@dataclass
class DryRunResult:
    """Result of a dry-run calculation.

    Attributes:
        output_dir: Directory containing generated files
        incar_path: Path to generated INCAR file
        poscar_path: Path to generated POSCAR file
        kpoints_path: Path to generated KPOINTS file
        potcar_path: Path to generated POTCAR file (may be None)
        success: Whether dry-run completed successfully
        error: Error message if success=False
    """

    output_dir: Path
    incar_path: Path
    poscar_path: Path
    kpoints_path: Path
    potcar_path: Optional[Path]
    success: bool
    error: Optional[str] = None

    def print_incar(self) -> None:
        """Print generated INCAR to stdout."""
        if self.incar_path.exists():
            print(self.incar_path.read_text())
        else:
            print(f"INCAR not found at {self.incar_path}")

    def print_poscar(self) -> None:
        """Print generated POSCAR to stdout."""
        if self.poscar_path.exists():
            print(self.poscar_path.read_text())
        else:
            print(f"POSCAR not found at {self.poscar_path}")

    def print_kpoints(self) -> None:
        """Print generated KPOINTS to stdout."""
        if self.kpoints_path.exists():
            print(self.kpoints_path.read_text())
        else:
            print(f"KPOINTS not found at {self.kpoints_path}")

    def get_incar_dict(self) -> Dict[str, Any]:
        """Parse INCAR file and return as dictionary."""
        if not self.incar_path.exists():
            return {}

        result = {}
        for line in self.incar_path.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('!'):
                continue
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip().upper()
                value = value.strip()
                # Remove inline comments
                if '#' in value:
                    value = value.split('#')[0].strip()
                if '!' in value:
                    value = value.split('!')[0].strip()
                result[key] = value
        return result

    def list_files(self) -> None:
        """List all generated files with sizes."""
        if not self.output_dir.exists():
            print(f"Output directory not found: {self.output_dir}")
            return

        print(f"\nFiles in {self.output_dir}:")
        for f in sorted(self.output_dir.iterdir()):
            size = f.stat().st_size
            if size < 1024:
                size_str = f"{size} B"
            elif size < 1024 * 1024:
                size_str = f"{size / 1024:.1f} KB"
            else:
                size_str = f"{size / (1024 * 1024):.1f} MB"
            print(f"  {f.name:<20} {size_str:>10}")


def dry_run_vasp(
    structure,
    builder_inputs: Dict[str, Any],
    code_label: str,
    output_dir: Optional[Union[str, Path]] = None,
    potential_family: str = 'PBE',
    potential_mapping: Optional[Dict[str, str]] = None,
    kpoints_spacing: Optional[float] = None,
    clean_submit_test: bool = True,
) -> DryRunResult:
    """Generate VASP input files without running calculation.

    Uses AiiDA's dry_run mode to prepare all inputs (INCAR, POSCAR,
    KPOINTS, POTCAR) without submitting to the scheduler.

    Args:
        structure: AiiDA StructureData or pymatgen Structure
        builder_inputs: Standard PS-TEROS builder_inputs dict
        code_label: Code label (e.g., 'vasp@localhost')
        output_dir: Where to save files (temp dir if None)
        potential_family: POTCAR family name
        potential_mapping: Element to POTCAR symbol mapping
        kpoints_spacing: Override kpoints_spacing from builder_inputs
        clean_submit_test: Remove submit_test folder after copying

    Returns:
        DryRunResult with paths to generated files

    Example:
        >>> from pymatgen.core import Structure, Lattice
        >>> structure = Structure(Lattice.cubic(4.0), ['Si', 'Si'],
        ...                       [[0, 0, 0], [0.25, 0.25, 0.25]])
        >>> builder_inputs = {
        ...     'parameters': {'incar': {'encut': 520, 'ismear': 0}},
        ...     'options': {'resources': {'num_machines': 1}},
        ...     'kpoints_spacing': 0.03,
        ... }
        >>> result = dry_run_vasp(structure, builder_inputs, 'vasp@localhost')
        >>> result.print_incar()
    """
    from aiida import orm
    from aiida.engine import run
    from aiida.plugins import WorkflowFactory

    # Setup output directory
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        output_dir = Path(tempfile.mkdtemp(prefix='vasp_dry_run_'))

    try:
        # Get VaspWorkChain
        VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
        builder = VaspWorkChain.get_builder()

        # Convert structure if needed
        if not isinstance(structure, orm.StructureData):
            structure = orm.StructureData(pymatgen=structure)

        # Set inputs
        builder.structure = structure
        builder.code = orm.load_code(code_label)

        # Parameters
        params = builder_inputs.get('parameters', {})
        builder.parameters = orm.Dict(dict=params)

        # Options
        options = builder_inputs.get('options', {})
        builder.options = orm.Dict(dict=options)

        # Potentials
        builder.potential_family = potential_family
        if potential_mapping:
            builder.potential_mapping = orm.Dict(dict=potential_mapping)
        else:
            # Auto-generate mapping from structure
            elements = list(set(str(s.symbol) for s in structure.get_pymatgen().species))
            builder.potential_mapping = orm.Dict(
                dict={el: el for el in elements}
            )

        # K-points
        ksp = kpoints_spacing or builder_inputs.get('kpoints_spacing')
        if ksp:
            builder.kpoints_spacing = ksp

        # Enable dry run
        builder.metadata.dry_run = True
        builder.metadata.store_provenance = False

        # Run dry-run
        try:
            run(builder)
        except Exception as e:
            # Dry run may raise but still generate files
            error_msg = str(e)

        # Find and copy generated files
        # AiiDA puts them in a submit_test folder
        submit_test = Path.cwd() / 'submit_test'

        if submit_test.exists():
            # Find the latest dry-run folder
            folders = sorted(
                [f for f in submit_test.iterdir() if f.is_dir()],
                key=lambda p: p.stat().st_mtime,
            )
            if folders:
                src = folders[-1]
                # Copy all files to output_dir
                for f in src.iterdir():
                    if f.is_file():
                        shutil.copy2(f, output_dir)

                # Clean up submit_test if requested
                if clean_submit_test:
                    shutil.rmtree(submit_test, ignore_errors=True)

        # Check for generated files
        incar_path = output_dir / 'INCAR'
        poscar_path = output_dir / 'POSCAR'
        kpoints_path = output_dir / 'KPOINTS'
        potcar_path = output_dir / 'POTCAR'

        success = incar_path.exists() or poscar_path.exists()

        return DryRunResult(
            output_dir=output_dir,
            incar_path=incar_path,
            poscar_path=poscar_path,
            kpoints_path=kpoints_path,
            potcar_path=potcar_path if potcar_path.exists() else None,
            success=success,
            error=None if success else "No input files generated",
        )

    except Exception as e:
        return DryRunResult(
            output_dir=output_dir,
            incar_path=output_dir / 'INCAR',
            poscar_path=output_dir / 'POSCAR',
            kpoints_path=output_dir / 'KPOINTS',
            potcar_path=None,
            success=False,
            error=str(e),
        )


def generate_incar_from_dict(
    incar_dict: Dict[str, Any],
    output_path: Optional[Union[str, Path]] = None,
) -> str:
    """Generate INCAR file content from dictionary.

    Does not require AiiDA - pure Python INCAR generation.

    Args:
        incar_dict: Dictionary of INCAR parameters
        output_path: Optional path to write INCAR file

    Returns:
        INCAR file content as string

    Example:
        >>> content = generate_incar_from_dict({
        ...     'ENCUT': 520,
        ...     'ISMEAR': 0,
        ...     'SIGMA': 0.05,
        ... })
        >>> print(content)
        ENCUT = 520
        ISMEAR = 0
        SIGMA = 0.05
    """
    lines = []

    # Sort keys for consistent output
    for key in sorted(incar_dict.keys()):
        value = incar_dict[key]

        # Format value based on type
        if isinstance(value, bool):
            value_str = '.TRUE.' if value else '.FALSE.'
        elif isinstance(value, (list, tuple)):
            value_str = ' '.join(str(v) for v in value)
        else:
            value_str = str(value)

        lines.append(f"{key.upper()} = {value_str}")

    content = '\n'.join(lines) + '\n'

    if output_path:
        Path(output_path).write_text(content)

    return content


def generate_kpoints_from_mesh(
    mesh: tuple,
    shift: tuple = (0, 0, 0),
    output_path: Optional[Union[str, Path]] = None,
) -> str:
    """Generate KPOINTS file content from mesh.

    Args:
        mesh: K-points mesh (kx, ky, kz)
        shift: K-points shift (default: Gamma-centered)
        output_path: Optional path to write KPOINTS file

    Returns:
        KPOINTS file content as string

    Example:
        >>> content = generate_kpoints_from_mesh((4, 4, 4))
        >>> print(content)
        Automatic mesh
        0
        Gamma
        4 4 4
        0 0 0
    """
    content = f"""Automatic mesh
0
Gamma
{mesh[0]} {mesh[1]} {mesh[2]}
{shift[0]} {shift[1]} {shift[2]}
"""

    if output_path:
        Path(output_path).write_text(content)

    return content
