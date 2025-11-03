#!/bin/bash
# Diagnostic script to find spack on cluster02

echo "=== Checking spack location on cluster02 ==="
echo ""

echo "1. Checking .bashrc for spack references:"
grep -n "spack" ~/.bashrc 2>&1 || echo "No spack found in ~/.bashrc"
echo ""

echo "2. Checking .bash_profile for spack references:"
grep -n "spack" ~/.bash_profile 2>&1 || echo "No .bash_profile or no spack references"
echo ""

echo "3. Searching for spack installation (may take a moment):"
find ~ -maxdepth 4 -type f -name "setup-env.sh" -path "*/spack/*" 2>/dev/null | head -5
echo ""

echo "4. Checking common spack locations:"
for path in ~/spack /opt/spack /usr/local/spack ~/.local/spack; do
    if [ -f "$path/share/spack/setup-env.sh" ]; then
        echo "FOUND: $path/share/spack/setup-env.sh"
    fi
done
echo ""

echo "5. Checking if spack is in PATH:"
which spack 2>&1 || echo "spack not in PATH"
echo ""

echo "6. Current directory:"
pwd
