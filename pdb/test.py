import os
import sys
import warnings
import numpy as np
import pandas as pd
from Bio.PDB import *
from Bio.PDB.MMCIF2Dict import *
from Bio.PDB.PDBExceptions import PDBConstructionWarning #ignorar warning (PDBConstructionWarning: WARNING: Chain B is discontinuous at line numeroDaLinha.)

warnings.simplefilter('ignore', PDBConstructionWarning)

os.chdir("/home/douglas_lima/pdb/testesCif")
file = open("/home/douglas_lima/pdb/testesCif/2wmg.cif")
file = open("/home/douglas_lima/pdb/testesCif/4of3.cif")

pdbParser = MMCIFParser()

structure = pdbParser.get_structure(file.name, file)
residueList = Selection.unfold_entities(structure, 'R')


mmcif_dict = MMCIF2Dict("/home/douglas_lima/pdb/testesCif/4of3.cif")

#oligo
entity_ids = mmcif_dict["_entity.id"]
entity_types = mmcif_dict["_entity.type"]
entity_descriptions = mmcif_dict["_entity.pdbx_description"]
entity_number_of_molecules = mmcif_dict["_entity.pdbx_number_of_molecules"]
entity_formula_weight = mmcif_dict["_entity.formula_weight"]

branch_entity_id = mmcif_dict["_pdbx_entity_branch.entity_id"]
branch_entity_type = mmcif_dict["_pdbx_entity_branch.type"]

#/////////////////////////////////////////////////////////////////////////////
#
# Se o branched retornar algo q nao for oligossacaridio da pra fazer assim:
#
#  dic = dict(zip(branch_entity_id, branch_entity_type))
# print([k for k, v in dic.items() if v == 'oligosaccharide'])

#print(mmcif_dict.keys())
#
#/////////////////////////////////////////////////////////////////////////////

# for i in range(len(entity_ids)) :
#     print(i)

#creates a dataframe
#Pandas package
#
residue_dict = {"ID": entity_ids, "type": entity_types, "description": entity_descriptions,"mol_num": entity_number_of_molecules, "formula_weight": entity_formula_weight, "entry ID": "4OF3"}
df = pd.DataFrame(data = residue_dict)
df = df[df.type == "branched"]


branch_list_entity_id = mmcif_dict["_pdbx_entity_branch_list.entity_id"]
branch_list_comp_id = mmcif_dict["_pdbx_entity_branch_list.comp_id"]
branch_list_comp_num = mmcif_dict["_pdbx_entity_branch_list.num"]

monossacaride_dict = {"comp ID": branch_list_comp_id, "entry ID": "4OF3", "oligossacaride": True, "entity_id": branch_list_entity_id, "comp_num": branch_list_comp_num, "mol_num": None}
monossacaride_df = pd.DataFrame(data = monossacaride_dict)
print(monossacaride_df)


#/////////////////////////////////////////////////////////////////////////////
#
# os.chdir("/home/douglas_lima/pdb/dicts")
# comp_dict_file = open("/home/douglas_lima/pdb/dicts/components.cif")

# print(comp_dict_file)


# comp_dict = MMCIF2Dict(comp_dict_file.name)
# comp_dict_file.close()
#
#/////////////////////////////////////////////////////////////////////////////

carbo_dict = pd.read_csv("/home/douglas_lima/pdb/dicts/CCD_carbohydrate_list.tsv", sep = "\t", header = None, names = ['carbo_id', 'REF'])

#print('SUM' in carbo_dict['carbo_id'].values)

#/////////////////////////////////////////////////////////////////////////////
#
# Verifica que o camp _chem_comp.id representa todos os componentes da estrutura
#
# entity_poly_seq_mon_ids = mmcif_dict["_entity_poly_seq.mon_id"]
# #print(entity_poly_seq_mon_ids)
# print(len(entity_poly_seq_mon_ids))
# print('GLY' in entity_poly_seq_mon_ids)

# chem_comp_ids = mmcif_dict["_chem_comp.id"]
# print(chem_comp_ids)
# i = 0
# for x in entity_poly_seq_mon_ids:
#     if(x in chem_comp_ids):
#         i = i + 1
# print(i)

#/////////////////////////////////////////////////////////////////////////////
#
# Teste da função set pra intersecção de duas listas
#
# a = [1, 2, 3, 4]
# b = [7, 6, 5, 4]
# a_set = set(a)
# b_set = set(b)
# print(set(a) & set(b))
# print(a_set & b_set)

# if(set(a) & set(b)):
#     print("Vrau")
#/////////////////////////////////////////////////////////////////////////////

#/////////////////////////////////////////////////////////////////////////////
#
# Criação da função pra filtrar estruturas com carboidratos
#
# def filter_containCarbo(fileNames):
    
#     filteredFileNames = []
#     carbo_dict = pd.read_csv("/home/douglas_lima/pdb/dicts/CCD_carbohydrate_list.tsv", sep = "\t", header = None, names = ['carbo_id', 'release_status']) # Release status values: https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp.pdbx_release_status.html

#     for fileName in fileNames:
        
#         mmcif_dict = MMCIF2Dict(fileName)
#         chem_components = mmcif_dict["_chem_comp.id"]

#         if(set(carbo_dict["carbo_id"].values) & set(chem_components)):
            
#             filteredFileNames.append(fileName)

#     return filteredFileNames

# fileNames = os.listdir("/home/douglas_lima/pdb/testesCif")
# print(filter_containCarbo(fileNames))
#
#/////////////////////////////////////////////////////////////////////////////