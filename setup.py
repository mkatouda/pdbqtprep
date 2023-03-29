from setuptools import setup

setup(
    name="pdbqtprep",
    version="0.1.1",
    install_requires=[
        'pyyaml', 'numpy', 'scipy', 'pdbfixer', 'openmm >= 7.1', 'meeko'
    ],
    entry_points={
        'console_scripts': [
            'pdbqtprep=pdbqtprep.pdbqtprep:main',
        ],
    },
    author="Michio Katouda",
    author_email="katouda@rist.or.jp",
    description="Autodock Vina protein and ligand PDBQT file prepare tool from PDB file",
    url="https://github.com/mkatouda/pdbqtprep",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires='>=3.7',
)
