#!/usr/bin/env python3
"""
Setup script for Teros package
"""

from setuptools import setup, find_packages

# Read version from __init__.py
with open("__init__.py", "r") as f:
    for line in f:
        if line.startswith("__version__"):
            version = line.split("=")[1].strip().strip('"').strip("'")
            break
    else:
        version = "0.1.0"

# Read description from __init__.py docstring
with open("__init__.py", "r") as f:
    content = f.read()
    start = content.find('"""') + 3
    end = content.find('"""', start)
    long_description = content[start:end].strip()

setup(
    name="teros",
    version=version,
    description="Ab Initio Atomistic Thermodynamics for Surface Energetics",
    long_description=long_description,
    long_description_content_type="text/plain",
    author="Teros Development Team",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "aiida-core",
        "aiida-workgraph", 
        "numpy",
        "matplotlib",
        "ase",
        "pymatgen",
        "rich",
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
)
