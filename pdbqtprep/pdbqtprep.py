#!/usr/bin/env python

import os
import sys
import argparse
import subprocess
import urllib
import pprint

import yaml
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from MolKit import Read
from AutoDockTools.MoleculePreparation import AD4ReceptorPreparation


def get_parser():
    class customHelpFormatter(argparse.ArgumentDefaultsHelpFormatter,
                              argparse.RawTextHelpFormatter):
        pass

    parser = argparse.ArgumentParser(
        formatter_class=customHelpFormatter,
        description="gromax protein-ligand MD trajectory analysis tools"
    )
    parser.add_argument(
        '-i', '--inp', type=str,
        help = 'yaml style input file, overwriting argument values'
    )
    parser.add_argument(
        '-p', '--pdbfile', type=str, default=None, dest='pdbfile',
        help = 'input PDB file'
    )
    parser.add_argument(
        '--pdbid', type=str, default=None, dest='pdbid',
        help='PDB id to retrieve from RCSB'
    )
    parser.add_argument(
        '--url', type=str, default=None, dest='url',
        help='URL to retrieve PDB from'
    )
    parser.add_argument(
        '-o', '--outbasename', type=str, default='output', dest='outbasename',
        help = 'basename of output PDB and PDBQT files'
    )
    parser.add_argument(
        '--add-atoms', default='all', dest='add-atoms', choices=('all', 'heavy', 'hydrogen', 'none'),
        help='which missing atoms to add: all, heavy, hydrogen, or none'
    )
    parser.add_argument(
        '--keep-heterogens', default='all', dest='keep-heterogens', choices=('all', 'water', 'none'),
        help='which heterogens to keep: all, water, or none'
    )
    parser.add_argument(
        '--replace-nonstandard', action='store_true', default=False, dest='replace-nonstandard',
        help='replace nonstandard residues with standard equivalents'
    )
    parser.add_argument(
        '--add-residues', action='store_true', default=False, dest='add-residues',
        help='add missing residues'
    )
    parser.add_argument(
        '--remove-chains', default=None, dest='remove-chains',
        help='specify list of indices or IDs of removing chains'
    )
    parser.add_argument(
        '--ph', type=float, default=7.4, dest='ph',
        help='the pH to use for adding missing hydrogens'
    )
    parser.add_argument(
        '--seed', type=int, default=None, dest='seed',
        help='Random seed'
    )
    parser.add_argument(
        '-v', '--verbose', action='store_true', dest='verbose',
        help='Print verbose output'
    )
    args = parser.parse_args()

    print(args)

    return args 

def set_config(args):
    # Read config yaml file
    if args.inp is not None and os.path.isfile(args.inp):
        with open(args.inp, 'r') as f:
            conf = yaml.safe_load(f)
    else:
        conf = {}

    # Set up default config values from program arguments
    conf_def = vars(args).copy()
    [conf.setdefault(k, v) for k, v in conf_def.items()]

    return conf

def pdbqtprep_main(conf):

    pdbid = conf['pdbid']
    pdburl = conf['url']
    pdbin_path = conf['pdbfile']
    outbasename = conf['outbasename']
    pdbout_path = os.path.abspath(outbasename + '_fixed.pdb')
    pdbqtout_path = os.path.abspath(outbasename + '_fixed.pdbqt')
    add_atoms = conf['add-atoms']
    keep_heterogens = conf['keep-heterogens']
    add_residues = conf['add-residues']
    ignore_terminal_missing_residues = conf['ignore_terminal_missing_residues']
    ph = conf['ph']
    removeChains = conf['remove-chains']
    seed = conf['seed']
    verbose = conf['verbose']

    if verbose:
        print('pdbin_path:', pdbin_path)
        print('pdbout_path:', pdbout_path)
        print('pdbqtout_path:', pdbqtout_path)

    if pdbid != None:
        print('Retrieving PDB "' + pdbid + '" from RCSB...')
        fixer = PDBFixer(pdbid=pdbid)
        pdburl = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdbid
        data = urllib.request.urlopen(pdburl).read()
        pdbin_path = os.path.abspath(outbasename + '.pdb')
        with open(pdbin_path, mode='wb') as f:
            f.write(data)
    elif pdburl != None:
        print('Retrieving PDB from URL "' + pdburl + '"...')
        fixer = PDBFixer(url=pdburl)
        data = urllib.request.urlopen(pdburl).read()
        pdbin_path = os.path.abspath(outbasename + '.pdb')
        with open(pdbin_path, mode='wb') as f:
            f.write(data)
    else:
        print('Retrieving PDB from file "' + pdbin_path + '"...')
        fixer = PDBFixer(filename=pdbin_path)

    if keep_heterogens == 'none':
        fixer.removeHeterogens(False)
    elif keep_heterogens == 'water':
        fixer.removeHeterogens(True)

    if isinstance(removeChains, list):
        print('removeChains:', removeChains)
        if isinstance(removeChains[0], int):
            fixer.removeChains(chainIndices=removeChains)
        elif isinstance(removeChains[0], str):
            fixer.removeChains(chainIds=removeChains)

    #with open(pdbin_path, 'w') as f:
    #    if fixer.source is not None:
    #        f.write("REMARK   1 PDBFIXER FROM: %s\n" % fixer.source)
    #    PDBFile.writeFile(fixer.topology, fixer.positions, f, True)

    if add_residues:
        print('Finding missing residues...')
        fixer.findMissingResidues()
        # if missing terminal residues shall be ignored, remove them from the dictionary
        if ignore_terminal_missing_residues:
            chains = list(fixer.topology.chains())
            keys = fixer.missingResidues.keys()
            for key in list(keys):
                chain = chains[key[0]]
                if key[1] == 0 or key[1] == len(list(chain.residues())):
                    del fixer.missingResidues[key]
        print('found ', len(fixer.missingResidues), 'missing regions')
        pprint.pprint(fixer.missingResidues)
    else:
        fixer.missingResidues = {}        

    print('Finding nonstandard residues...')
    fixer.findNonstandardResidues()
    print('Replacing nonstandard residues...')
    fixer.replaceNonstandardResidues()

    fixer.findMissingAtoms()
    if add_atoms not in ('all', 'heavy'):
        fixer.missingAtoms = {}
        fixer.missingTerminals = {}
    print('found ', len(fixer.missingAtoms), 'missing atoms')
    pprint.pprint(fixer.missingAtoms)
    print('found ', len(fixer.missingTerminals), 'missing terminals')
    pprint.pprint(fixer.missingTerminals)
    fixer.addMissingAtoms(seed=seed)

    if add_atoms in ('all', 'hydrogen'):
        print('Adding missing hydrogens...')
        fixer.addMissingHydrogens(ph)

    with open(pdbout_path, 'w') as f:
        if fixer.source is not None:
            f.write("REMARK   1 PDBFIXER FROM: %s\n" % fixer.source)
        PDBFile.writeFile(fixer.topology, fixer.positions, f, True)

    mode = 'automatic'
    repairs = ''
    charges_to_add = 'gasteiger'
    cleanup  = "nphs_lps_waters_nonstdres"
    preserved = {}
    delete_single_nonstd_residues = None
    dictionary = None

    #mols = Read(pdbin_path)
    mols = Read(pdbout_path)
    mol = mols[0]

    unique_atom_names = False
    if unique_atom_names:  # added to simplify setting up covalent dockings 8/2014
        for at in mol.allAtoms:
            if mol.allAtoms.get(at.name) >1:
                at.name = at.name + str(at._uniqIndex +1)
        if verbose:
            print("renamed %d atoms: each newname is the original name of the atom plus its (1-based) uniqIndex" %(len(mol.allAtoms)))

    mol.buildBondsByDistance()
    alt_loc_ats = mol.allAtoms.get(lambda x: "@" in x.name)
    len_alt_loc_ats = len(alt_loc_ats)
    if len_alt_loc_ats:
        print("WARNING!", mol.name, "has",len_alt_loc_ats, ' alternate location atoms!\nUse prepare_pdb_split_alt_confs.py to create pdb files containing a single conformation.\n')

    RPO = AD4ReceptorPreparation(mol, mode, repairs, charges_to_add, 
                        cleanup, outputfilename=pdbqtout_path,
                        preserved=preserved, 
                        delete_single_nonstd_residues=delete_single_nonstd_residues,
                        dict=dictionary)    

def main():
    args = get_parser()
    if args.verbose: print(args)

    conf = set_config(args)

    print('======= Input configulations =======')
    for k, v in conf.items():
        print('{}: {}'.format(k, v))
    print('====================================')

    pdbqtprep_main(conf)

if __name__ == '__main__':
    main()
