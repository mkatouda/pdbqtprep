# pdbqtprep

Integrated program of PDBFixer and AutoDockTools to prepare protein and ligand model of    
AutoDock Vina docking and/or molecular dynamics simulation simulation.

## Licence

This package is distributed under the MIT License.

## Required softwares

1. python: 3.7 or later
2. pyyaml (https://pyyaml.org/)
3. pdbfixer (https://github.com/openmm/pdbfixer)
4. rdkit (https://github.com/rdkit/rdkit)
5. scipy (https://github.com/scipy/scipy)
5. meeko (https://github.com/forlilab/Meeko)
6. autodocktools_py3 (https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3)

## Installation

## Installation (Required)

- Create conda virtual environment  
```
conda create -n py38-pdbqtprep python=3.8  
conda activate py38-pdbqtprep  
```

- Run the following command to install required conda packages  
```
conda install -c conda-forge pyyaml pdbfixer rdkit scipy  
```

- Install meeko from PyPI
```
pip install meeko
```

- Install autodocktools_py3 from github
```
pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3.git
```

- Install pdbqtprep from github  
```
pip install git+https://github.com/mkatouda/pdbqtprep.git
```

- Install pdbqtprep from local repository  
```
git clone https://github.com/mkatouda/pdbqtprep.git
cd pdbqreprep
pip install .
```

## Command usage

```
usage: pdbqtprep [-h] [-i INP] [-p PDBFILE] [--pdbid PDBID] [--url URL] [-o OUTBASENAME] [--add-atoms {all,heavy,hydrogen,none}] [--keep-heterogens {all,water,none}]
                 [--replace-nonstandard] [--add-residues] [--remove-chains REMOVE-CHAINS] [--ph PH] [--seed SEED] [--ligid LIGID] [-v]

protein and ligand coordinate files from a complex pdb file

options:
  -h, --help            show this help message and exit
  -i INP, --inp INP     yaml style input file, overwriting argument values (default: None)
  -p PDBFILE, --pdbfile PDBFILE
                        input PDB file (default: None)
  --pdbid PDBID         PDB id to retrieve from RCSB (default: None)
  --url URL             URL to retrieve PDB from (default: None)
  -o OUTBASENAME, --outbasename OUTBASENAME
                        basename of output PDB and PDBQT files (default: output)
  --add-atoms {all,heavy,hydrogen,none}
                        which missing atoms to add: all, heavy, hydrogen, or none (default: all)
  --keep-heterogens {all,water,none}
                        which heterogens to keep: all, water, or none (default: all)
  --replace-nonstandard
                        replace nonstandard residues with standard equivalents (default: False)
  --add-residues        add missing residues (default: False)
  --remove-chains REMOVE-CHAINS
                        specify list of indices or IDs of removing chains (default: None)
  --ph PH               the pH to use for adding missing hydrogens (default: 7.4)
  --seed SEED           Random seed (default: None)
  --ligid LIGID         Ligand ID in PDB file (default: None)
  -v, --verbose         Print verbose output (default: False)
```

## Exmaples of command line usage

```
pdbqtprep --pdbid 3POZ --ligid 03P -o 3POZ --add-atoms all --keep-heterogens none --replace-nonstandard --add-residues --ph 7.4
```

## Exmaples of yaml input usage

Prepare input yaml file input.yml:

```
pdbid: '3POZ'
ligid: '03P'
outbasename: '3POZ'
keep-heterogens: 'none'
add-atoms: 'all'
replace-nonstandard: True
add-residues: True
ph: 7.4
ignore_terminal_missing_residues: True
verbose: True
```

Then, run pdbqtprep in command line:

```
pdbqtprep -i input.yml
```

You can get protein models (3POZ_pro.pdb, 3POZ_pro.pdbqt) and ligand models (3POZ_lig.mol, 3POZ_lig.pdbqt).    

Keywards of yaml file are the same in the name of command line options.  
See above explation of command line options.  

## Author

Michio Katouda (katouda@rist.or.jp)  
