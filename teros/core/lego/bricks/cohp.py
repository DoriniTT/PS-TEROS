"""COHP brick for the lego module.

Handles COHP (Crystal Orbital Hamilton Population) analysis using LOBSTER.
"""

import os
import glob
import shutil
import tempfile
import subprocess
from pathlib import Path

from aiida import orm
from aiida.common.links import LinkType
from aiida_workgraph import task
from .connections import COHP_PORTS as PORTS  # noqa: F401


def validate_stage(stage: dict, stage_names: set) -> None:
    """Validate a COHP stage configuration.

    Args:
        stage: Stage configuration dict.
        stage_names: Set of stage names defined so far (before this stage).

    Raises:
        ValueError: If validation fails.
    """
    name = stage['name']

    if 'cohp_from' not in stage:
        raise ValueError(
            f"Stage '{name}': cohp stages require 'cohp_from' "
            f"(name of stage with VASP outputs)"
        )

    cohp_from = stage['cohp_from']
    if cohp_from not in stage_names:
        raise ValueError(
            f"Stage '{name}' cohp_from='{cohp_from}' must reference "
            f"a previous stage name"
        )


def create_stage_tasks(wg, stage, stage_name, context):
    """Create COHP stage tasks in the WorkGraph.

    Args:
        wg: WorkGraph to add tasks to.
        stage: Stage configuration dict.
        stage_name: Unique stage identifier.
        context: Dict with shared context.

    Returns:
        Dict with task references for later stages.
    """
    stage_tasks = context['stage_tasks']
    stage_types = context['stage_types']
    stages = context['stages']

    cohp_from = stage['cohp_from']

    # Get the remote folder and retrieved FolderData from the referenced stage
    ref_stage_type = stage_types.get(cohp_from, 'vasp')
    if ref_stage_type == 'vasp':
        remote_socket = stage_tasks[cohp_from]['vasp'].outputs.remote_folder
        retrieved_socket = stage_tasks[cohp_from]['vasp'].outputs.retrieved
    else:
        raise ValueError(
            f"COHP stage '{stage_name}' cohp_from='{cohp_from}' must "
            f"reference a VASP stage (got type='{ref_stage_type}')"
        )

    # Resolve structure: prefer output structure (from relaxation),
    # fall back to input structure (for SCF with NSW=0)
    cohp_from_stage = stage_tasks[cohp_from]
    cohp_from_incar = next(
        (s.get('incar', {}) for s in stages if s['name'] == cohp_from),
        {},
    )
    if cohp_from_incar.get('nsw', 0) > 0:
        stage_structure = cohp_from_stage['vasp'].outputs.structure
    else:
        stage_structure = cohp_from_stage.get(
            'input_structure',
            cohp_from_stage['vasp'].outputs.structure,
        )

    # Prepare LOBSTER settings
    lobsterin_settings = {
        'cohp_start_energy': stage.get('cohp_start_energy', -15.0),
        'cohp_end_energy': stage.get('cohp_end_energy', 10.0),
        'basis_set': stage.get('basis_set', 'pbeVaspFit2015'),
        'bond_cutoff': stage.get('bond_cutoff', 3.0),
        'save_projection': stage.get('save_projection', True),
        'compute_coop': stage.get('compute_coop', True),
    }

    # Add LOBSTER analysis task
    cohp_task = wg.add_task(
        run_lobster_analysis,
        name=f'cohp_{stage_name}',
        remote_folder=remote_socket,
        retrieved=retrieved_socket,
        structure=stage_structure,
        lobsterin_settings=orm.Dict(dict=lobsterin_settings),
    )

    return {
        'cohp': cohp_task,
        'structure': stage_structure,
    }


def expose_stage_outputs(wg, stage_name, stage_tasks_result, namespace_map=None):
    """Expose COHP stage outputs on the WorkGraph.

    Args:
        wg: WorkGraph instance.
        stage_name: Unique stage identifier.
        stage_tasks_result: Dict returned by create_stage_tasks.
        namespace_map: Dict mapping output group to namespace string,
                      e.g. {'main': 'stage1'}. If None, falls back to
                      flat naming with stage_name prefix.
    """
    cohp_task = stage_tasks_result['cohp']

    if namespace_map is not None:
        ns = namespace_map['main']
        setattr(wg.outputs, f'{ns}.cohp.cohp_data', cohp_task.outputs.cohp_data)
        setattr(wg.outputs, f'{ns}.cohp.icohp', cohp_task.outputs.icohp)
        setattr(wg.outputs, f'{ns}.cohp.cohpcar', cohp_task.outputs.cohpcar)
        setattr(wg.outputs, f'{ns}.cohp.doscar', cohp_task.outputs.doscar)
    else:
        setattr(wg.outputs, f'{stage_name}_cohp_data', cohp_task.outputs.cohp_data)
        setattr(wg.outputs, f'{stage_name}_icohp', cohp_task.outputs.icohp)
        setattr(wg.outputs, f'{stage_name}_cohpcar', cohp_task.outputs.cohpcar)
        setattr(wg.outputs, f'{stage_name}_doscar', cohp_task.outputs.doscar)


def get_stage_results(wg_node, wg_pk: int, stage_name: str,
                      namespace_map: dict = None) -> dict:
    """Extract results from a COHP stage in a sequential workflow.

    Args:
        wg_node: The WorkGraph ProcessNode.
        wg_pk: WorkGraph PK.
        stage_name: Name of the COHP stage.
        namespace_map: Dict mapping output group to namespace string,
                      e.g. {'main': 'stage1'}. If None, uses flat naming.

    Returns:
        Dict with keys: cohp_data, lobster_files, pk, stage, type.
    """
    result = {
        'cohp_data': None,
        'lobster_files': {},
        'pk': wg_pk,
        'stage': stage_name,
        'type': 'cohp',
    }

    if hasattr(wg_node, 'outputs'):
        outputs = wg_node.outputs

        if namespace_map is not None:
            ns = namespace_map['main']
            stage_ns = getattr(outputs, ns, None)
            brick_ns = getattr(stage_ns, 'cohp', None) if stage_ns is not None else None
            if brick_ns is not None:
                if hasattr(brick_ns, 'cohp_data'):
                    cohp_data_node = brick_ns.cohp_data
                    if hasattr(cohp_data_node, 'get_dict'):
                        result['cohp_data'] = cohp_data_node.get_dict()
                for file_key in ('icohp', 'cohpcar', 'doscar'):
                    if hasattr(brick_ns, file_key):
                        result['lobster_files'][file_key] = getattr(brick_ns, file_key)
        else:
            # Flat naming fallback
            # COHP data Dict
            cohp_data_attr = f'{stage_name}_cohp_data'
            if hasattr(outputs, cohp_data_attr):
                cohp_data_node = getattr(outputs, cohp_data_attr)
                if hasattr(cohp_data_node, 'get_dict'):
                    result['cohp_data'] = cohp_data_node.get_dict()

            # LOBSTER output files
            for file_key in ('icohp', 'cohpcar', 'doscar'):
                file_attr = f'{stage_name}_{file_key}'
                if hasattr(outputs, file_attr):
                    result['lobster_files'][file_key] = getattr(outputs, file_attr)

    # Fallback: traverse links
    if result['cohp_data'] is None:
        _extract_cohp_stage_from_workgraph(wg_node, stage_name, result)

    return result


def _extract_cohp_stage_from_workgraph(
    wg_node, stage_name: str, result: dict
) -> None:
    """Extract COHP stage results by traversing WorkGraph links.

    Args:
        wg_node: The WorkGraph ProcessNode.
        stage_name: Name of the COHP stage.
        result: Result dict to populate (modified in place).
    """
    if not hasattr(wg_node, 'base'):
        return

    cohp_task_name = f'cohp_{stage_name}'

    called_calc = wg_node.base.links.get_outgoing(link_type=LinkType.CALL_CALC)
    for link in called_calc.all():
        child_node = link.node
        link_label = link.link_label

        if cohp_task_name in link_label or link_label == cohp_task_name:
            created = child_node.base.links.get_outgoing(link_type=LinkType.CREATE)
            for out_link in created.all():
                out_label = out_link.link_label
                out_node = out_link.node

                if out_label == 'cohp_data' and hasattr(out_node, 'get_dict'):
                    result['cohp_data'] = out_node.get_dict()
                elif out_label in ('icohp', 'cohpcar', 'doscar'):
                    result['lobster_files'][out_label] = out_node
            break


def print_stage_results(index: int, stage_name: str, stage_result: dict) -> None:
    """Print formatted results for a COHP stage.

    Args:
        index: 1-based stage index.
        stage_name: Name of the stage.
        stage_result: Result dict from get_stage_results.
    """
    print(f"  [{index}] {stage_name} (COHP)")

    if stage_result['cohp_data'] is not None:
        cohp_data = stage_result['cohp_data']
        bonds = cohp_data.get('bonds', [])
        print(f"      Bonds analyzed: {len(bonds)}")

        # Print summary of strongest bonds (top 5 by |ICOHP|)
        if bonds:
            sorted_bonds = sorted(
                bonds,
                key=lambda b: abs(b.get('icohp', 0.0)),
                reverse=True
            )
            print("      Strongest bonding interactions (top 5 by |ICOHP|):")
            for i, bond in enumerate(sorted_bonds[:5], 1):
                atom1 = bond.get('atom1', '?')
                atom2 = bond.get('atom2', '?')
                distance = bond.get('distance', 0.0)
                icohp = bond.get('icohp', 0.0)
                print(
                    f"        {i}. {atom1}-{atom2}: "
                    f"d={distance:.3f} Å, ICOHP={icohp:.4f} eV"
                )

    lobster_files = stage_result.get('lobster_files', {})
    if lobster_files:
        file_names = ', '.join(
            f"{k} (PK {v.pk})" for k, v in lobster_files.items()
        )
        print(f"      LOBSTER files: {file_names}")


# ─── LOBSTER calcfunction tasks ──────────────────────────────────────────────


@task.calcfunction(outputs=['cohp_data', 'icohp', 'cohpcar', 'doscar'])
def run_lobster_analysis(
    remote_folder: orm.RemoteData,
    retrieved: orm.FolderData,
    structure: orm.StructureData,
    lobsterin_settings: orm.Dict,
) -> dict:
    """
    Run LOBSTER COHP analysis on VASP outputs.

    This calcfunction:
    1. Extracts required VASP files from retrieved FolderData
    2. Copies WAVECAR from remote folder (if available locally)
    3. Generates lobsterin configuration file
    4. Runs the lobster binary
    5. Parses ICOHPLIST.lobster and returns outputs

    Args:
        remote_folder: RemoteData from VASP calculation with WAVECAR
        retrieved: FolderData from VASP with CONTCAR, POTCAR, INCAR, DOSCAR
        structure: StructureData for the calculation
        lobsterin_settings: Dict with LOBSTER configuration parameters

    Returns:
        dict with:
            - 'cohp_data': orm.Dict with parsed ICOHPLIST.lobster data
            - 'icohp': orm.SinglefileData for ICOHPLIST.lobster
            - 'cohpcar': orm.SinglefileData for COHPCAR.lobster
            - 'doscar': orm.SinglefileData for DOSCAR.lobster
    """
    # Find the lobster binary
    lobster_path = Path.home() / '.local' / 'bin' / 'lobster'
    if not lobster_path.exists():
        # Try system PATH
        lobster_path = shutil.which('lobster')
        if lobster_path is None:
            raise FileNotFoundError(
                "LOBSTER binary not found. Install it to ~/.local/bin/lobster or add to PATH"
            )
        lobster_path = Path(lobster_path)

    # Create a temporary working directory
    tmpdir = tempfile.mkdtemp(prefix='lobster_')
    try:
        # Extract required files from retrieved FolderData
        required_retrieved = ['CONTCAR', 'POTCAR', 'INCAR', 'DOSCAR']
        for fname in required_retrieved:
            try:
                content = retrieved.get_object_content(fname, mode='rb')
                filepath = os.path.join(tmpdir, fname)
                with open(filepath, 'wb') as f:
                    f.write(content)
            except (FileNotFoundError, OSError) as e:
                raise FileNotFoundError(
                    f"File '{fname}' not found in retrieved FolderData (PK {retrieved.pk}). "
                    f"Ensure the SCF stage has these files in the retrieve list. Error: {e}"
                )

        # Try to get WAVECAR from retrieved first, fall back to remote
        wavecar_path = os.path.join(tmpdir, 'WAVECAR')
        try:
            wavecar_content = retrieved.get_object_content('WAVECAR', mode='rb')
            with open(wavecar_path, 'wb') as f:
                f.write(wavecar_content)
        except (FileNotFoundError, OSError):
            # Try to copy from remote folder (this only works if running on same machine)
            # In production, WAVECAR should be in retrieved list
            if hasattr(remote_folder, 'get_remote_path'):
                remote_wavecar = os.path.join(remote_folder.get_remote_path(), 'WAVECAR')
                if os.path.exists(remote_wavecar):
                    shutil.copy(remote_wavecar, wavecar_path)
                else:
                    raise FileNotFoundError(
                        f"WAVECAR not found in retrieved (PK {retrieved.pk}) "
                        f"or remote folder (PK {remote_folder.pk}). "
                        f"Add WAVECAR to the VASP stage retrieve list."
                    )
            else:
                raise FileNotFoundError(
                    f"WAVECAR not found in retrieved FolderData (PK {retrieved.pk}). "
                    f"Add WAVECAR to the VASP stage retrieve list."
                )

        # Create lobsterin file
        settings = lobsterin_settings.get_dict()
        lobsterin_content = _generate_lobsterin(settings)
        lobsterin_path = os.path.join(tmpdir, 'lobsterin')
        with open(lobsterin_path, 'w') as f:
            f.write(lobsterin_content)

        # Run lobster binary
        lobster_result = subprocess.run(
            [str(lobster_path)],
            cwd=tmpdir,
            capture_output=True,
            text=True,
        )

        if lobster_result.returncode != 0:
            raise RuntimeError(
                f"LOBSTER analysis failed (return code {lobster_result.returncode}).\n"
                f"stdout: {lobster_result.stdout}\n"
                f"stderr: {lobster_result.stderr}"
            )

        # Parse ICOHPLIST.lobster
        icohplist_path = os.path.join(tmpdir, 'ICOHPLIST.lobster')
        if not os.path.exists(icohplist_path):
            raise FileNotFoundError(
                f"ICOHPLIST.lobster not produced by LOBSTER. "
                f"Check LOBSTER output:\n{lobster_result.stdout}"
            )

        cohp_data = _parse_icohplist(icohplist_path)

        # Build output dict
        outputs = {}
        outputs['cohp_data'] = orm.Dict(dict=cohp_data)

        # Store LOBSTER output files as SinglefileData
        lobster_files = {
            'icohp': 'ICOHPLIST.lobster',
            'cohpcar': 'COHPCAR.lobster',
            'doscar': 'DOSCAR.lobster',
        }

        for key, filename in lobster_files.items():
            filepath = os.path.join(tmpdir, filename)
            if os.path.exists(filepath):
                outputs[key] = orm.SinglefileData(file=filepath)
            else:
                # Create empty placeholder if file not generated
                with open(filepath, 'w') as f:
                    f.write(f"# {filename} not generated by LOBSTER\n")
                outputs[key] = orm.SinglefileData(file=filepath)

        return outputs

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def _generate_lobsterin(settings: dict) -> str:
    """
    Generate lobsterin configuration file content.

    Args:
        settings: Dict with LOBSTER parameters

    Returns:
        String content for lobsterin file
    """
    lines = []

    # Energy range for COHP
    cohp_start = settings.get('cohp_start_energy', -15.0)
    cohp_end = settings.get('cohp_end_energy', 10.0)
    lines.append(f'COHPstartEnergy  {cohp_start:.1f}')
    lines.append(f'COHPendEnergy     {cohp_end:.1f}')

    # Basis set
    basis_set = settings.get('basis_set', 'pbeVaspFit2015')
    lines.append(f'basisSet          {basis_set}')

    # COHP generator (automatic bond detection)
    bond_cutoff = settings.get('bond_cutoff', 3.0)
    lines.append(f'cohpGenerator     from 0.1 to {bond_cutoff:.1f} orbitalwise')

    # Optional settings
    if settings.get('save_projection', True):
        lines.append('saveProjectionToFile')

    if settings.get('compute_coop', True):
        lines.append('createFatband COOPCAR.lobster')

    # Add custom lines if provided
    custom_lines = settings.get('custom_lines', [])
    if custom_lines:
        lines.extend(custom_lines)

    return '\n'.join(lines) + '\n'


def _parse_icohplist(filepath: str) -> dict:
    """
    Parse ICOHPLIST.lobster file.

    Args:
        filepath: Path to ICOHPLIST.lobster file

    Returns:
        dict with bonds list and summary info
    """
    bonds = []

    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Find the header line
    header_idx = None
    for i, line in enumerate(lines):
        if 'No.' in line and 'Distance' in line and 'ICOHP' in line:
            header_idx = i
            break

    if header_idx is None:
        return {'bonds': bonds, 'total_bonds': 0}

    # Parse bond entries
    for line in lines[header_idx + 1:]:
        line = line.strip()
        if not line or line.startswith('-'):
            continue

        # Format: No.  Atom1  Atom2  Length  ICOHP  IntegrationFermi
        parts = line.split()
        if len(parts) < 6:
            continue

        try:
            bond = {
                'number': int(parts[0]),
                'atom1': parts[1],
                'atom2': parts[2],
                'distance': float(parts[3]),
                'icohp': float(parts[4]),
                'integration_fermi': float(parts[5]),
            }
            bonds.append(bond)
        except (ValueError, IndexError):
            continue

    return {
        'bonds': bonds,
        'total_bonds': len(bonds),
    }
