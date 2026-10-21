"""Vibrational surface energy of Cu(100) from phonons on a q-point mesh, without AiiDA (ASE's EMT gives the forces).

The pure-Python layer of the phonopy mode (``psteros.phonons``): ``phonopy_displacements`` gives phonopy's
symmetry-reduced displacements of a structure in its supercell, any calculator gives their forces (here EMT; in the
graphs a DFT static or the fine-tuned model of the MLIP route), and ``phonopy_modes`` the modes on a q-point mesh. The
slab vibrates as a whole, in an in-plane supercell with a 2D mesh, and

    gamma_vib(T) = [F_vib(slab) - n F_vib(bulk)] / 2A

for a slab of n bulk cells. The script converges it with the slab thickness and the q-point mesh.

    python cu100_emt.py          # needs phonopy, ase and pymatgen; about 30 s
"""

from __future__ import annotations

from ase.build import bulk, fcc100
from ase.calculators.emt import EMT
from ase.optimize import BFGS
from pymatgen.io.ase import AseAtomsAdaptor

import psteros

TEMPERATURES = (300.0, 600.0, 900.0)


def forces(structures):
    """EMT forces (eV/A) of pymatgen structures."""

    result = []
    for structure in structures:
        atoms = AseAtomsAdaptor.get_atoms(structure)
        atoms.calc = EMT()
        result.append(atoms.get_forces())
    return result


def modes(label, structure, config):
    """phonopy displacements -> forces -> modes on the mesh of ``config`` for ``label``."""

    supercell, mesh = config.supercell(label, structure), config.mesh(label, structure)
    displacements = psteros.phonopy_displacements(structure, supercell)
    perfect, displaced = displacements.structures(structure)
    result = psteros.phonopy_modes(structure, displacements, forces(displaced),
                                   reference_forces_ev_per_angstrom=forces([perfect])[0], mesh=mesh)
    print(f"  {label:8s} supercell {supercell}, {len(displacements)} displacements, mesh {mesh}: "
          f"{result.stability.summary()}")
    return result.vibrations


def relaxed_slab(layers):
    atoms = fcc100("Cu", size=(1, 1, layers), a=3.59, vacuum=8.0, periodic=True)
    atoms.calc = EMT()
    BFGS(atoms, logfile=None).run(fmax=1e-4)
    return AseAtomsAdaptor.get_structure(atoms)


def main() -> None:
    cu = AseAtomsAdaptor.get_structure(bulk("Cu", "fcc", a=3.59))
    for mesh_length in (20.0, 40.0):
        config = psteros.PhonopyConfig(supercell_length_angstrom=10.0, mesh_length_angstrom=mesh_length, slabs=("slab",))
        print(f"mesh length {mesh_length:g} A")
        bulk_modes = modes("bulk", cu, config)
        for layers in (3, 5, 7):
            slab = relaxed_slab(layers)
            slab_modes = modes("slab", slab, config)
            area = slab.lattice.a * slab.lattice.b
            gamma = [1000 * (slab_modes.free_energy_ev(t) - layers * bulk_modes.free_energy_ev(t)) / (2 * area)
                     for t in TEMPERATURES]
            print(f"    {layers} layers: gamma_vib " + ", ".join(f"{value:6.2f}" for value in gamma)
                  + f" meV/A^2 at {', '.join(f'{t:g}' for t in TEMPERATURES)} K")


if __name__ == "__main__":
    main()
