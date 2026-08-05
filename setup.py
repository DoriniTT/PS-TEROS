from setuptools import setup, find_packages

setup(
    name="psteros",
    version="1.0.0",
    description="Reproducible surface thermodynamics with AiiDA",
    author="Thiago T. Dorini",
    author_email="your.email@example.com",
    packages=find_packages(),
    install_requires=[
        "aiida-core",
        "aiida-workgraph",
        "aiida-quantumespresso",
        "aiida-pseudo",
        "pymatgen",
        "numpy",
    ],
    extras_require={
        "vasp": [
            "aiida-vasp>=5,<6",
        ],
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
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.11",
)
