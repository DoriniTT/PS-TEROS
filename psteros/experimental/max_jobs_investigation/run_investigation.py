#!/usr/bin/env python3
"""
Run all max_number_jobs investigation tests.

This script runs all three approaches and documents the findings.
"""

import sys
import os

# Add current directory to path for imports
sys.path.insert(0, os.path.dirname(__file__))

from aiida import load_profile
load_profile()

print("""
================================================================================
MAX_NUMBER_JOBS INVESTIGATION
================================================================================

Problem: max_number_jobs on main WorkGraph doesn't limit concurrent VASP
calculations in nested sub-workgraphs created by @task.graph decorators.

We're testing three approaches to solve this:
  1. Parameter Passing - Add max_number_jobs as function parameter
  2. Decorator Configuration - Pass max_number_jobs to @task.graph decorator
  3. Framework Properties - Check if max_number_jobs cascades from parent

================================================================================
""")

results = {}

# Test Approach 1: Parameter Passing
print("\n" + "="*80)
print("APPROACH 1: Parameter Passing")
print("="*80)

try:
    from test_approach_1_parameter import test_parameter_passing
    pk1 = test_parameter_passing()
    results['approach_1'] = {
        'status': 'completed',
        'pk': pk1,
        'description': 'Tested adding max_number_jobs as function parameter'
    }
except Exception as e:
    print(f"\n[ERROR] Approach 1 failed: {e}")
    import traceback
    traceback.print_exc()
    results['approach_1'] = {
        'status': 'failed',
        'error': str(e)
    }

# Test Approach 2: Decorator Configuration
print("\n" + "="*80)
print("APPROACH 2: Decorator Configuration")
print("="*80)

try:
    from test_approach_2_decorator import test_decorator_configuration
    pk2 = test_decorator_configuration()
    results['approach_2'] = {
        'status': 'completed' if pk2 else 'not_supported',
        'pk': pk2,
        'description': 'Tested passing max_number_jobs to @task.graph() decorator'
    }
except Exception as e:
    print(f"\n[ERROR] Approach 2 failed: {e}")
    import traceback
    traceback.print_exc()
    results['approach_2'] = {
        'status': 'failed',
        'error': str(e)
    }

# Test Approach 3: Framework Properties
print("\n" + "="*80)
print("APPROACH 3: Framework Properties")
print("="*80)

try:
    from test_approach_3_framework import test_framework_cascade, test_framework_properties
    test_framework_properties()  # Inspection first
    pk3 = test_framework_cascade()
    results['approach_3'] = {
        'status': 'completed',
        'pk': pk3,
        'description': 'Tested if max_number_jobs cascades from parent to nested workgraphs'
    }
except Exception as e:
    print(f"\n[ERROR] Approach 3 failed: {e}")
    import traceback
    traceback.print_exc()
    results['approach_3'] = {
        'status': 'failed',
        'error': str(e)
    }

# Print summary
print("\n" + "="*80)
print("INVESTIGATION SUMMARY")
print("="*80)

for approach, result in results.items():
    print(f"\n{approach.upper()}:")
    print(f"  Status: {result['status']}")
    if 'pk' in result and result['pk']:
        print(f"  WorkGraph PK: {result['pk']}")
        print(f"  Check with: verdi process show {result['pk']}")
    if 'description' in result:
        print(f"  Description: {result['description']}")
    if 'error' in result:
        print(f"  Error: {result['error']}")

print("\n" + "="*80)
print("NEXT STEPS")
print("="*80)
print("""
1. Check the process reports for timing information:
   - Do calculations run 2 at a time (respects max_number_jobs)?
   - Or do all 5 run simultaneously (ignores max_number_jobs)?

2. Inspect the workgraph structure:
   - verdi process show <PK>
   - Look for nested workgraphs
   - Check their max_number_jobs settings

3. Document which approach (if any) works

Run: verdi process list -a -p 1
""")
