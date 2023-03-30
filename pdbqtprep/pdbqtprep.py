#!/usr/bin/env python

import os
import sys
import argparse
import subprocess
import urllib
import requests
import json
import pprint

import yaml
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from MolKit import Read
from AutoDockTools.MoleculePreparation import AD4ReceptorPreparation
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from meeko import MoleculePreparation


def get_parser():
    class customHelpFormatter(argparse.ArgumentDefaultsHelpFormatter,
                              argparse.RawTextHelpFormatter):
        pass

    parser = argparse.ArgumentParser(
        formatter_class=customHelpFormatter,
        description="protein and ligand coordinate files from a complex pdb file"
    )
    parser.add_argument(
        '-i', '--inp', type=str,
        help = 'yaml style input file, overwriting argument values'
    )
    parser.add_argument(
        '-p', '--pdbfile', type=str, default=None, dest='pdbfile',
        help = 'input protein PDB file'
    )
    parser.add_argument(
        '-l', '--ligfile', type=str, default=None, dest='ligfile',
        help = 'input ligand file (mol, sdf)'
    )
    parser.add_argument(
        '--pdbid', type=str, default=None, dest='pdbid',
        help='PDB ID to retrieve from RCSB'
    )
    parser.add_argument(
        '--ligid', type=str, default=None, dest='ligid',
        help='Ligand residue name to retrieve from RCSB'
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
    return parser.parse_args()

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

def pdbqt_pro(pdbid, pdburl, pdbin_path, pro_pdbout_path, pro_pdbqtout_path,
              add_atoms='all', keep_heterogens='all', add_residues=False,
              ignore_terminal_missing_residues=True,
              removeChains=None, ph=7.4, seed=None, verbose=False):

    if pdbid != None:
        print('Retrieving PDB "' + pdbid + '" from RCSB...')
        fixer = PDBFixer(pdbid=pdbid)
        pdburl = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdbid
        data = urllib.request.urlopen(pdburl).read()
        pdbin_path = os.path.abspath(pdbid + '.pdb')
        with open(pdbin_path, mode='wb') as f:
            f.write(data)
    elif pdburl != None:
        print('Retrieving PDB from URL "' + pdburl + '"...')
        fixer = PDBFixer(url=pdburl)
        data = urllib.request.urlopen(pdburl).read()
        pdbin_path = os.path.abspath(pdburl + '.pdb')
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

    with open(pro_pdbout_path, 'w') as f:
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
    mols = Read(pro_pdbout_path)
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
        print("WARNING!", mol.name, "has", len_alt_loc_ats, ' alternate location atoms!\nUse prepare_pdb_split_alt_confs.py to create pdb files containing a single conformation.\n')

    RPO = AD4ReceptorPreparation(
        mol, mode, repairs, charges_to_add, 
        cleanup,
        outputfilename=pro_pdbqtout_path,
        preserved=preserved, 
        delete_single_nonstd_residues=delete_single_nonstd_residues,
        dict=dictionary
    )

    print('final pdbin_path:', pdbin_path)
    return pdbin_path

def get_lig_from_PDB(entry_id, ligid, encoding='sdf', target_auth_asym_id='A',
                     verbose=False):
    url = f"https://data.rcsb.org/rest/v1/core/entry/{entry_id}"
    risposta = requests.get(url)
    parsed = json.loads(risposta.content.decode())
    parsed_lig = parsed["rcsb_entry_container_identifiers"]["non_polymer_entity_ids"]
    #print(parsed_lig)

    try:
        parsed_binding_affinity = parsed["rcsb_binding_affinity"]
        #print(parsed_binding_affinity, len(parsed_binding_affinity))
        #print(parsed_binding_affinity[0])
        binding_comp_id = parsed_binding_affinity[0]["comp_id"]
    except:
        if ligid is not None:
            binding_comp_id = ligid
        else:
            binding_comp_id = None
    print('binding_comp_id:', binding_comp_id)

    parsed_lig_dict = {}
    for lig in parsed_lig:
        url = f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{entry_id}/{lig}"
        risposta = requests.get(url)
        parsed = json.loads(risposta.content.decode())
        parsed_lig_dict[lig] = [
            parsed["pdbx_entity_nonpoly"]["comp_id"],
            dict.fromkeys(parsed["rcsb_nonpolymer_entity_container_identifiers"]["asym_ids"]),
        ]
    if verbose: print(parsed_lig_dict)

    for lig in parsed_lig_dict:
        for i, chain in enumerate(parsed_lig_dict[lig][1]):
            url = f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity_instance/{entry_id}/{chain}"
            risposta = requests.get(url)
            parsed = json.loads(risposta.content.decode())
            parsed_lig_dict[lig][1][chain] = [
                parsed["rcsb_nonpolymer_entity_instance_container_identifiers"]["auth_asym_id"],
                parsed["rcsb_nonpolymer_entity_instance_container_identifiers"]["auth_seq_id"]
            ]
    if verbose: print(parsed_lig_dict)

    ligout_path = None
    for lig in parsed_lig_dict:
        comp_id = parsed_lig_dict[lig][0]
        for chain in parsed_lig_dict[lig][1]:
            auth_asym_id = parsed_lig_dict[lig][1][chain][0]
            seq_id = parsed_lig_dict[lig][1][chain][1]
            url = f"https://models.rcsb.org/v1/{entry_id}/ligand?auth_seq_id={seq_id}&label_asym_id={chain}&encoding={encoding}"
            risposta = requests.get(url, allow_redirects=True)
            out_path = f"{entry_id}_lig_{comp_id}_{chain}.{encoding}"
            #out_path = f"{entry_id}_{chain}_{comp_id}.{encoding}"
            print(lig, comp_id, seq_id, chain, auth_asym_id, url, out_path)
            with open(out_path, 'wb') as f:
                f.write(risposta.content)

            if auth_asym_id == target_auth_asym_id and comp_id == binding_comp_id:
                ligout_path = out_path
                lig_molout_path = f"{entry_id}_lig_{comp_id}_{chain}_addH.mol"
                lig_pdbqtout_path = f"{entry_id}_lig_{comp_id}_{chain}_addH.pdbqt"
                gen_lig_pdbqt(out_path, lig_molout_path, lig_pdbqtout_path, verbose=verbose)
                #ligout_path = f"{entry_id}_lig.{encoding}"
                #print(lig, comp_id, seq_id, chain, auth_asym_id, url, out_path)
                #with open(ligout_path, 'wb') as f:
                #    f.write(risposta.content)

    return ligout_path

def gen_lig_pdbqt(ligin_path, lig_molout_path, lig_pdbqtout_path, verbose=False):
    mol = Chem.MolFromMolFile(ligin_path)
    mol = Chem.RemoveHs(mol)
    mol = Chem.AddHs(mol, addCoords=True)
    Chem.MolToMolFile(mol, lig_molout_path)
    mol_conf = mol.GetConformer(-1)
    com = list(rdMolTransforms.ComputeCentroid(mol_conf))

    print('COM of ligand:', com)
    if verbose: print(Chem.MolToMolBlock(mol))

    config = MoleculePreparation.get_defaults_dict()
    if verbose: print('MoleculePreparation config:', config)
    preparator = MoleculePreparation.from_config(config)
    preparator.prepare(mol)
    mol_pdbqt = preparator.write_pdbqt_string()
    with open(lig_pdbqtout_path, 'w') as f:
        f.write(mol_pdbqt)

def pdbqt_lig(pdbid, ligid, ligin_path, lig_molout_path, lig_pdbqtout_path, verbose=False):

    if pdbid != None:
        ligin_path = get_lig_from_PDB(
            pdbid, ligid, encoding='mol',
            target_auth_asym_id='A', 
            verbose=verbose
        )
    print('ligin_path:', ligin_path)

    gen_lig_pdbqt(ligin_path, lig_molout_path, lig_pdbqtout_path, verbose=verbose)

def pdbqtprep_main(conf):

    pdbid = conf['pdbid']
    pdburl = conf['url']
    pdbin_path = conf['pdbfile']
    ligid = conf['ligid']
    ligin_path = conf['ligfile']
    outbasename = conf['outbasename']
    add_atoms = conf['add-atoms']
    keep_heterogens = conf['keep-heterogens']
    add_residues = conf['add-residues']
    ignore_terminal_missing_residues = conf['ignore_terminal_missing_residues']
    removeChains = conf['remove-chains']
    ph = conf['ph']
    seed = conf['seed']
    verbose = conf['verbose']

    pro_pdbout_path = os.path.abspath(outbasename + '_pro.pdb')
    pro_pdbqtout_path = os.path.abspath(outbasename + '_pro.pdbqt')
    lig_molout_path = os.path.abspath(outbasename + '_lig.mol')
    lig_pdbout_path = os.path.abspath(outbasename + '_lig.pdb')
    lig_pdbqtout_path = os.path.abspath(outbasename + '_lig.pdbqt')

    print('pro_pdbout_path:', pro_pdbout_path)
    print('pro_pdbqtout_path:', pro_pdbqtout_path)

    pdbin_path = pdbqt_pro(pdbid, pdburl, pdbin_path, pro_pdbout_path, pro_pdbqtout_path,
                           add_atoms=add_atoms, keep_heterogens=keep_heterogens,
                           add_residues=add_residues,
                           ignore_terminal_missing_residues=ignore_terminal_missing_residues,
                           removeChains=removeChains, ph=ph, seed=seed, verbose=verbose)

    print('lig_molout_path:', lig_molout_path)
    print('lig_pdbqtout_path:', lig_pdbqtout_path)

    pdbqt_lig(pdbid, ligid, ligin_path, lig_molout_path, lig_pdbqtout_path, verbose=verbose)
                  
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
