from setuptools import setup, find_packages

setup(
    name="teros",
    version="0.1.0",
    description="Ab Initio Atomistic Thermodynamics for Surface Energetics",
    author="Thiago T. Dorini",
    author_email="your.email@example.com",
    packages=find_packages(),
    install_requires=[
        "aiida-core",
        "aiida-workgraph",
        "pymatgen",
        "numpy",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
)