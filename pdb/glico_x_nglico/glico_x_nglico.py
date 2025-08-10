from Bio import PDB
import os
import pandas as pd
from pathlib import Path
from tqdm import tqdm

DATA_DIR_PATH = Path('/home/douglas/carboanalysis/data/unzipped/')
CARBO_LIST_PATH = Path('/home/douglas/carboanalysis/carboanalysis/pdb/dicts/CCD_carbohydrate_list.tsv')

def get_all_covalent_contacts(mmcif_file):
    parser = PDB.MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", DATA_DIR_PATH / mmcif_file)

    carbo_dict = pd.read_csv(CARBO_LIST_PATH, sep="\t", header=None, names=['carbo_id', 'release_status'])
    carbohydrate_residues = carbo_dict["carbo_id"].values

    protein_atoms = []
    carbohydrate_atoms = []
    carbohydrate_chains = set()

    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()
                if resname in carbohydrate_residues:
                    carbohydrate_chains.add(chain)
                    for atom in residue.get_atoms():
                        carbohydrate_atoms.append((residue, atom))
                elif PDB.is_aa(residue):
                    for atom in residue.get_atoms():
                        protein_atoms.append((residue, atom))

    results = set()
    covalent_chains = set()
    for protein_residue, protein_atom in protein_atoms:
        for carb_residue, carb_atom in carbohydrate_atoms:
            dist = protein_atom - carb_atom
            # (entry, carb_resname, carb_residue_number, carb_chain_id)
            results.add((mmcif_file.replace('.cif', ''),  # sem extensão
                            carb_residue.get_resname(),
                            carb_residue.get_id()[1],
                            carb_residue.get_parent().get_id()))
            if dist <= 2.0:
                covalent_chains.add(carb_residue.get_parent().get_id())

    if results:
        carbohydrate_residues_df = pd.DataFrame(results, columns=['entry', 'carb_resname', 'carb_residue_number', 'carb_chain_id'])
        carbohydrate_residues_df['covalent'] = False
        carbohydrate_residues_df.loc[carbohydrate_residues_df['carb_chain_id'].isin(covalent_chains), 'covalent'] = True
        
        return carbohydrate_residues_df
    else:
        return pd.DataFrame(columns=['entry', 'carb_resname', 'carb_residue_number', 'carb_chain_id', 'covalent'])


entry_df = pd.read_csv('/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/carbo_entrys_res_owab_filtered_nd.txt', names=['entry'])
entry_list = entry_df['entry'].values

all_results = []
for entry in tqdm(entry_list, desc='Caracterizando Interações'):
    df = get_all_covalent_contacts(entry)
    all_results.append(df)

final_df = pd.concat(all_results, ignore_index=True)
final_df.to_csv('carbohydrate_interactions.csv', index=False)