from setuptools import setup, find_packages

setup(
    name="teros",
    version="0.2.0",
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
    extras_require={
        "docs": [
            "sphinx>=4.0.0",
            "sphinx-rtd-theme>=1.0.0",
            "recommonmark",
            "sphinxcontrib-bibtex",
        ],
        "dev": [
            "pytest",
            "pytest-cov",
            "black",
            "flake8",
        ]
    },
    entry_points={
        'aiida.calculations': [
            'teros.aimd_vasp = teros.core.lego.calcs.aimd_vasp:AimdVaspCalculation',
        ],
        'aiida.workflows': [
            'teros.aimd_vasp_workchain = teros.core.lego.calcs.aimd_vasp:AimdVaspWorkChain',
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.11",
)