#!/usr/bin/env python
"""
Monitor relaxation energy calculation and display results.

Usage:
    python monitor_relaxation_energy.py <PK>
"""

import sys
import time
from aiida import load_profile, orm

def main():
    if len(sys.argv) < 2:
        print("Usage: python monitor_relaxation_energy.py <PK>")
        sys.exit(1)
    
    pk = int(sys.argv[1])
    
    load_profile()
    
    print(f"\n{'='*80}")
    print(f"MONITORING WORKGRAPH PK: {pk}")
    print(f"{'='*80}\n")
    
    node = orm.load_node(pk)
    
    print(f"WorkGraph: {node.label}")
    print(f"State: {node.process_state.value}")
    print(f"\nWaiting for completion...")
    print(f"(Press Ctrl+C to stop monitoring, workflow will continue running)\n")
    
    last_state = None
    
    try:
        while not node.is_terminated:
            current_state = node.process_state.value
            
            if current_state != last_state:
                print(f"[{time.strftime('%H:%M:%S')}] State: {current_state}")
                last_state = current_state
            
            time.sleep(10)
            node = orm.load_node(pk)  # Reload to get latest state
        
        print(f"\n[{time.strftime('%H:%M:%S')}] State: {node.process_state.value}")
        print(f"\n{'='*80}")
        print(f"WORKFLOW COMPLETED!")
        print(f"{'='*80}\n")
        
        if node.process_state.value == 'finished':
            print("✓ Workflow finished successfully!")
            
            # Check if relaxation_energies output exists
            if hasattr(node.outputs, 'relaxation_energies'):
                print(f"\n{'-'*80}")
                print("RELAXATION ENERGIES")
                print(f"{'-'*80}")
                
                for label, energy in node.outputs.relaxation_energies.items():
                    print(f"  {label}: {energy.value:+.4f} eV")
                
                print(f"\n{'-'*80}")
                print("ALL OUTPUTS AVAILABLE")
                print(f"{'-'*80}")
                print(f"  • slab_structures:         {len(node.outputs.slab_structures)} slabs")
                print(f"  • unrelaxed_slab_energies: {len(node.outputs.unrelaxed_slab_energies)} energies")
                print(f"  • relaxed_slabs:           {len(node.outputs.relaxed_slabs)} structures")
                print(f"  • slab_energies:           {len(node.outputs.slab_energies)} energies")
                print(f"  • relaxation_energies:     {len(node.outputs.relaxation_energies)} energies")
            else:
                print("\n! relaxation_energies output not found")
                print(f"Available outputs: {list(node.outputs)}")
        else:
            print(f"✗ Workflow ended with state: {node.process_state.value}")
            print(f"\nCheck report:")
            print(f"  verdi process report {pk}")
        
    except KeyboardInterrupt:
        print(f"\n\nMonitoring stopped (workflow continues running)")
        print(f"\nCheck status:")
        print(f"  verdi process show {pk}")
        print(f"\nResume monitoring:")
        print(f"  python monitor_relaxation_energy.py {pk}")

if __name__ == '__main__':
    main()
