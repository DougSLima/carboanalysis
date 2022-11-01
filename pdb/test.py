from ctypes import sizeof
from itertools import count
import os
import sys
import warnings
import numpy as np
import pandas as pd
from Bio.PDB import *
from Bio.PDB.MMCIF2Dict import *
from Bio.PDB.PDBExceptions import PDBConstructionWarning #ignorar warning (PDBConstructionWarning: WARNING: Chain B is discontinuous at line numeroDaLinha.)
from collections import defaultdict

warnings.simplefilter('ignore', PDBConstructionWarning)

os.chdir("/home/douglas_lima/pdb/testesCif")
file = open("/home/douglas_lima/pdb/testesCif/2wmg.cif")
file = open("/home/douglas_lima/pdb/testesCif/4of3.cif")

pdbParser = MMCIFParser()

structure = pdbParser.get_structure(file.name, file)
residueList = Selection.unfold_entities(structure, 'R')

#/////////////////////////////////////////////////////////////////////////////
#
# Pega os resíduos com hetero flag difernte de " "(Amino acidos) e "W" (Moléculas de água)
#
# for residue in residueList:
#     if residue.get_id()[0] != " " and residue.get_id()[0] != "W":
#         print(residue)
#/////////////////////////////////////////////////////////////////////////////

#/////////////////////////////////////////////////////////////////////////////
#
mmcif_dict = MMCIF2Dict(file.name)

entity_ids = mmcif_dict["_entity.id"]
entity_types = mmcif_dict["_entity.type"]
entity_descriptions = mmcif_dict["_entity.pdbx_description"]

#print(mmcif_dict.keys())
#
#/////////////////////////////////////////////////////////////////////////////

# for i in range(len(entity_ids)) :
#     print(i)

#creates a dataframe
#Pandas package
#
residue_dict = {"ID": entity_ids, "type": entity_types, "description": entity_descriptions, "Entry ID": "4OF3"}
df = pd.DataFrame(data = residue_dict)
#print(df)


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

entity_poly_seq_mon_ids = mmcif_dict["_entity_poly_seq.mon_id"]
print(entity_poly_seq_mon_ids)
print(len(entity_poly_seq_mon_ids))
print('GLY' in entity_poly_seq_mon_ids)

chem_comp_ids = mmcif_dict["_chem_comp.id"]









