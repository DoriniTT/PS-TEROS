#!/home/thiagotd/envs/psteros/bin/python
"""
Workaround script to access DOS and Bands results when BandsWorkChain fails.

The BandsWorkChain may fail to expose outputs correctly due to an aiida-vasp bug,
but the actual calculations complete successfully. This script shows how to 
access the results directly from the child workchains.

Usage:
    python workaround_access_results.py <workgraph_pk>
"""

import sys
from aiida import load_profile, orm

def access_bands_results_workaround(workgraph_pk):
    """
    Access bands and DOS results from a WorkGraph that used BandsWorkChain.
    
    Even if BandsWorkChain fails, we can access results from its child workchains.
    """
    load_profile('psteros')
    
    # Load the main workgraph
    wg = orm.load_node(workgraph_pk)
    print(f"WorkGraph: {wg.label} (PK: {wg.pk})")
    print(f"State: {wg.process_state}")
    print()
    
    # Find the BandsWorkChain
    bands_wc = None
    for called in wg.base.links.get_outgoing().all():
        node = called.node
        if 'BandsWorkChain' in str(type(node).__name__):
            bands_wc = node
            break
    
    if not bands_wc:
        print("❌ Could not find BandsWorkChain")
        return None
    
    print(f"BandsWorkChain: {bands_wc.label} (PK: {bands_wc.pk})")
    print(f"State: {bands_wc.process_state}")
    print()
    
    # Get the child workchains
    scf_wc = None
    bands_calc_wc = None
    dos_calc_wc = None
    
    for link in bands_wc.base.links.get_outgoing(link_type='call_work').all():
        node = link.node
        label = link.link_label
        
        if label == 'scf':
            scf_wc = node
        elif label == 'bs':
            bands_calc_wc = node
        elif label == 'dos':
            dos_calc_wc = node
    
    print("=" * 80)
    print("RESULTS")
    print("=" * 80)
    
    results = {}
    
    # Access relaxed/primitive structure
    if 'primitive_structure' in bands_wc.outputs:
        structure = bands_wc.outputs.primitive_structure
        results['structure'] = structure
        print(f"\n✓ Primitive Structure: {structure.uuid}")
        print(f"  Formula: {structure.get_formula()}")
        print(f"  Number of atoms: {len(structure.sites)}")
    
    # Access seekpath parameters
    if 'seekpath_parameters' in bands_wc.outputs:
        seekpath = bands_wc.outputs.seekpath_parameters
        results['seekpath_parameters'] = seekpath
        print(f"\n✓ Seekpath Parameters: {seekpath.uuid}")
        seekpath_dict = seekpath.get_dict()
        if 'path' in seekpath_dict:
            print(f"  High-symmetry path: {' → '.join(seekpath_dict['path'])}")
    
    # Access SCF results
    if scf_wc:
        print(f"\n✓ SCF Calculation: {scf_wc.uuid} (PK: {scf_wc.pk})")
        if 'misc' in scf_wc.outputs:
            misc = scf_wc.outputs.misc.get_dict()
            if 'total_energies' in misc and 'energy_extrapolated' in misc['total_energies']:
                energy = misc['total_energies']['energy_extrapolated']
                results['scf_energy'] = energy
                print(f"  Energy: {energy:.6f} eV")
    
    # Access Band Structure results
    if bands_calc_wc:
        print(f"\n✓ Band Structure Calculation: {bands_calc_wc.uuid} (PK: {bands_calc_wc.pk})")
        if 'bands' in bands_calc_wc.outputs:
            bands = bands_calc_wc.outputs.bands
            results['bands'] = bands
            print(f"  Bands: {bands.uuid}")
            try:
                band_arrays = bands.get_bands()
                print(f"  Shape: {band_arrays.shape} (kpoints, bands)")
                print(f"  Number of k-points: {band_arrays.shape[0]}")
                print(f"  Number of bands: {band_arrays.shape[1]}")
            except Exception as e:
                print(f"  Could not get band shape: {e}")
    
    # Access DOS results
    if dos_calc_wc:
        print(f"\n✓ DOS Calculation: {dos_calc_wc.uuid} (PK: {dos_calc_wc.pk})")
        # DOS is stored as 'bands' output in VaspWorkChain
        if 'bands' in dos_calc_wc.outputs:
            dos_bands = dos_calc_wc.outputs.bands
            results['dos'] = dos_bands
            print(f"  DOS data: {dos_bands.uuid}")
            try:
                # For DOS, the BandsData contains energy vs DOS
                dos_arrays = dos_bands.get_bands()
                print(f"  Shape: {dos_arrays.shape}")
                # Try to get the array data
                arrays = dos_bands.get_array('kpoints')
                print(f"  K-points array shape: {arrays.shape}")
            except Exception as e:
                print(f"  DOS array info: {e}")
    
    print("\n" + "=" * 80)
    print("EXPORT INSTRUCTIONS")
    print("=" * 80)
    
    print("\nTo export the data in Python:")
    print(f"```python")
    print(f"from aiida.orm import load_node")
    print()
    
    if 'bands' in results:
        pk = results['bands'].pk
        print(f"# Load band structure")
        print(f"bands = load_node({pk})")
        print(f"band_arrays = bands.get_bands()  # Shape: (n_kpoints, n_bands)")
        print(f"kpoints = bands.get_kpoints()    # K-point coordinates")
        print()
    
    if 'dos' in results:
        pk = results['dos'].pk
        print(f"# Load DOS data")
        print(f"dos = load_node({pk})")
        print(f"dos_data = dos.get_bands()  # DOS values")
        print()
    
    if 'structure' in results:
        pk = results['structure'].pk
        print(f"# Load optimized structure")
        print(f"structure = load_node({pk})")
        print(f"atoms = structure.get_ase()  # Convert to ASE")
        print(f"atoms.write('optimized_structure.cif')")
        print()
    
    print(f"```")
    
    return results


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python workaround_access_results.py <workgraph_pk>")
        print("\nExample:")
        print("  python workaround_access_results.py 29346")
        sys.exit(1)
    
    try:
        wg_pk = int(sys.argv[1])
    except ValueError:
        print(f"Error: '{sys.argv[1]}' is not a valid PK")
        sys.exit(1)
    
    results = access_bands_results_workaround(wg_pk)
    
    if results:
        print("\n" + "=" * 80)
        print("✓ Successfully accessed all results!")
        print("=" * 80)
