# Phonons on a q-point mesh (phonopy)

`cu100_emt.py` uses the AiiDA-free layer of the phonopy mode of psteros (`psteros.phonons`, see
`docs/source/vibrations.rst`): phonopy's displacements of the bulk and of Cu(100) slabs, their forces from ASE's EMT
potential, the modes on a q-point mesh, and the vibrational surface energy
`gamma_vib(T) = [F_vib(slab) - n F_vib(bulk)] / 2A` for 3, 5 and 7 layers and two meshes. It needs phonopy
(`pip install "psteros[phonopy]"`), ase and pymatgen, and runs in about half a minute.

In the graphs the same mode is `VibrationsConfig(phonopy=PhonopyConfig(...))`: `build_vibrations_workgraph` (DFT),
`dft_vibrations` and the `vibrations` of the MLIP route (`examples/mlip_routes/sno2_charge_neutral_mlip.py --phonopy`).
