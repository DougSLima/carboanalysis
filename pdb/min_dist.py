from Bio import PDB
import os
import pandas as pd

def get_min_distance_mmcif(mmcif_file):
    parser = PDB.MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", mmcif_file)
    
    # Lista de resíduos comuns de carboidratos no PDB
    carbo_dict = pd.read_csv("/home/douglas/carboanalysis/carboanalysis/pdb/dicts/CCD_carbohydrate_list.tsv", sep = "\t", header = None, names = ['carbo_id', 'release_status'])
    carbohydrate_residues = carbo_dict["carbo_id"].values

    protein_atoms = []
    carbohydrate_atoms = []

    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()
                if resname in carbohydrate_residues:
                    carbohydrate_atoms.extend(residue.get_atoms())
                elif PDB.is_aa(residue):  # Verifica se é um aminoácido
                    protein_atoms.extend(residue.get_atoms())

    # Calcular a menor distância
    min_distance = float("inf")
    closest_atoms = None

    for protein_atom in protein_atoms:
        for carb_atom in carbohydrate_atoms:
            dist = protein_atom - carb_atom  # Distância entre os átomos
            if dist < min_distance:
                min_distance = dist
                closest_atoms = (protein_atom, carb_atom)

    if closest_atoms:
        print(f"Menor distância: {min_distance:.3f} Å")
        print(f"Entre os átomos: {closest_atoms[0]} e {closest_atoms[1]}")
        print(closest_atoms[0].get_parent())
        print(closest_atoms[0].get_parent().get_resname())
        print(closest_atoms[0].get_parent().get_id()[1])
        print(closest_atoms[1].get_parent())
        print(closest_atoms[1].get_parent().get_resname())
        print(closest_atoms[1].get_parent().get_id()[1])
    else:
        print("Nenhuma interação encontrada.")

# Exemplo de uso

os.chdir("/home/douglas/carboanalysis/data/unzipped")
get_min_distance_mmcif("2wmg.cif")
