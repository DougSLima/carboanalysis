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
#file = open("/home/douglas_lima/pdb/testesCif/1v6u.cif")
#file = open("/home/douglas_lima/pdb/testesCif/5ebw.cif")

pdbParser = MMCIFParser()

structure = pdbParser.get_structure(file.name, file)
residueList = Selection.unfold_entities(structure, 'R')


mmcif_dict = MMCIF2Dict(file.name)

#entity
entity_ids = mmcif_dict["_entity.id"]
entity_types = mmcif_dict["_entity.type"]
entity_descriptions = mmcif_dict["_entity.pdbx_description"]
entity_number_of_molecules = mmcif_dict["_entity.pdbx_number_of_molecules"]
entity_formula_weight = mmcif_dict["_entity.formula_weight"]

if "_pdbx_entity_branch.entity_id" in mmcif_dict:
    branch_entity_id = mmcif_dict["_pdbx_entity_branch.entity_id"]
    branch_entity_type = mmcif_dict["_pdbx_entity_branch.type"]

#/////////////////////////////////////////////////////////////////////////////
#
# Se o branched retornar algo q nao for oligossacaridio da pra fazer assim:
#
# dic = dict(zip(branch_entity_id, branch_entity_type))
# print([k for k, v in dic.items() if v == 'oligosaccharide'])

#print(mmcif_dict.keys())
#
#/////////////////////////////////////////////////////////////////////////////

# for i in range(len(entity_ids)) :
#     print(i)

#creates a dataframe
#Pandas package
#
residue_dict = {"id": entity_ids, "type": entity_types, "description": entity_descriptions,"mol_num": entity_number_of_molecules, "formula_weight": entity_formula_weight, "entry_id": mmcif_dict["_entry.id"][0]}
df = pd.DataFrame(data = residue_dict)
#print(df)
df_branched = df[df.type == "branched"]
#print(df_branched)

if "_pdbx_entity_branch_list.entity_id" in mmcif_dict:
    branch_list_entity_id = mmcif_dict["_pdbx_entity_branch_list.entity_id"]
    branch_list_comp_id = mmcif_dict["_pdbx_entity_branch_list.comp_id"]
    branch_list_comp_num = mmcif_dict["_pdbx_entity_branch_list.num"]

    mol_nums = []
    for entity_id in branch_list_entity_id:
        mol_nums.append(entity_number_of_molecules[int(entity_id) - 1]) # transforma o vetor de ids pra inteiro e retira 1 pq o id de entidades começa em 1 e o indice do vetor começa em 0

    monosaccharide_dict = {"comp_id": branch_list_comp_id, "entry_id": mmcif_dict["_entry.id"][0], "oligossacaride": True, "entity_id": branch_list_entity_id, "comp_num": branch_list_comp_num, "mol_num": mol_nums}
    monosaccharide_df = pd.DataFrame(data = monosaccharide_dict)

    #print(monosaccharide_df)

if "_pdbx_entity_nonpoly.entity_id" in mmcif_dict:
    nonpoly_entity_id = mmcif_dict["_pdbx_entity_nonpoly.entity_id"]
    nonpoly_entity_comp_id = mmcif_dict["_pdbx_entity_nonpoly.comp_id"]

    mol_nums = []
    for entity_id in nonpoly_entity_id:
        mol_nums.append(entity_number_of_molecules[int(entity_id) - 1]) # transforma o vetor de ids pra inteiro e retira 1 pq o id de entidades começa em 1 e o indice do vetor começa em 0

    monosaccharide_dict_2 = {"comp_id": nonpoly_entity_comp_id, "entry_id": mmcif_dict["_entry.id"][0], "oligossacaride": False, "entity_id":  nonpoly_entity_id, "comp_num": None, "mol_num": mol_nums}
    monosaccharide_df_2 = pd.DataFrame(data = monosaccharide_dict_2)

    carbo_dict = pd.read_csv("/home/douglas_lima/pdb/dicts/CCD_carbohydrate_list.tsv", sep = "\t", header = None, names = ['carbo_id', 'REF'])

    monosaccharide_df_2 = monosaccharide_df_2[monosaccharide_df_2.comp_id.isin(carbo_dict["carbo_id"].values)]
    #print(monosaccharide_df_2)

monosaccharide = pd.concat([monosaccharide_df, monosaccharide_df_2], ignore_index=True)
#print(monosaccharide)

#Chem comp identifier

identifier_comp_id = mmcif_dict["_pdbx_chem_comp_identifier.comp_id"]
identifier_type = mmcif_dict["_pdbx_chem_comp_identifier.type"]
identifier_identifier = mmcif_dict["_pdbx_chem_comp_identifier.identifier"]

identifier_dict = {"comp_id": identifier_comp_id, "type": identifier_type, "identifier": identifier_identifier}
identifier_df = pd.DataFrame(data=identifier_dict)
#print(identifier_df[(identifier_df.comp_id == 'AHR') & (identifier_df.type == 'IUPAC CARBOHYDRATE SYMBOL')]["identifier"])
#print(identifier_df)

commom_names = []
iupac_symbols = []

for comp_id in monosaccharide.comp_id:
    commom_names.append(identifier_df[(identifier_df.comp_id == comp_id) & (identifier_df.type == 'COMMON NAME')]["identifier"].values[0]) 
    iupac_symbols.append(identifier_df[(identifier_df.comp_id == comp_id) & (identifier_df.type == 'IUPAC CARBOHYDRATE SYMBOL')]["identifier"].values[0]) 

monosaccharide["commom_name"] = commom_names
monosaccharide["iupac_symbol"] = iupac_symbols

# x = monosaccharide[monosaccharide.comp_id.isin(['GLU', 'NAG'])]
# print(x)
#print(monosaccharide)

#monosaccharide.to_csv(path_or_buf="/home/douglas_lima/pdb/monosaccharide.csv")

def filter_containCarbo(fileNames):
    
    filteredFileNames = []
    filteredEntries = []
    carbo_dict = pd.read_csv("/home/douglas_lima/pdb/dicts/CCD_carbohydrate_list.tsv", sep = "\t", header = None, names = ['carbo_id', 'release_status']) # Release status values: https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp.pdbx_release_status.html

    for fileName in fileNames:
        
        mmcif_dict = MMCIF2Dict(fileName)
        chem_components = mmcif_dict["_chem_comp.id"]

        if(set(carbo_dict["carbo_id"].values) & set(chem_components)):
            
            filteredFileNames.append(fileName)
            filteredEntries.append(mmcif_dict["_entry.id"][0])

    #Escreve um arquivo .txt com o nome das entradas cujas estrutuaras possuam carbohidratos
    with open("/home/douglas_lima/pdb/dataframes/carbo_entrys.txt", "w") as file:
        for entry in filteredEntries:
            #escreve cada entry_id em uma nova linha
            file.write("%s\n" % entry)
            
    return filteredFileNames

#Separa os monossacarídeos e oligossacarídeos em dois pd.DataFrames e gera os respectivos arquivos .csv
def separate(fileNames):
    
    fileNames = filter_containCarbo(fileNames)

    oligosaccharide_df = pd.DataFrame()
    olig_monosaccharide_df = pd.DataFrame()
    monosaccharide_df = pd.DataFrame()
    monosaccharides = pd.DataFrame()

    for fileName in fileNames:

        #cria um dicionário a partir do arquivo .cif
        mmcif_dict = MMCIF2Dict(fileName)
        
        #entity data
        entity_ids = mmcif_dict["_entity.id"]
        entity_types = mmcif_dict["_entity.type"]
        entity_descriptions = mmcif_dict["_entity.pdbx_description"]
        entity_number_of_molecules = mmcif_dict["_entity.pdbx_number_of_molecules"]
        entity_formula_weight = mmcif_dict["_entity.formula_weight"]

        #oligosaccharide DataFrame
        if "_pdbx_entity_branch.entity_id" in mmcif_dict:

            entity_dict = {"id": entity_ids, "entry_id": mmcif_dict["_entry.id"][0], "type": entity_types, "description": entity_descriptions,"mol_num": entity_number_of_molecules, "formula_weight": entity_formula_weight}
            entity_df = pd.DataFrame(data = entity_dict)
            
            branch_entity_id = mmcif_dict["_pdbx_entity_branch.entity_id"]
            branch_entity_type = mmcif_dict["_pdbx_entity_branch.type"]

            #usa a lista de ids de entidade e de tipo pra verificar quais entidades são oligossacarídeos
            olig_dic = dict(zip(branch_entity_id, branch_entity_type))

            entity_df = entity_df[entity_df.id.isin([k for k, v in olig_dic.items() if v == 'oligosaccharide'])]
            
            oligosaccharide_df = pd.concat([oligosaccharide_df, entity_df], ignore_index=True)

        #monosaccharide from oligosaccharide DataFrame
        if "_pdbx_entity_branch_list.entity_id" in mmcif_dict:

            branch_list_entity_id = mmcif_dict["_pdbx_entity_branch_list.entity_id"]
            branch_list_comp_id = mmcif_dict["_pdbx_entity_branch_list.comp_id"]
            branch_list_comp_num = mmcif_dict["_pdbx_entity_branch_list.num"]

            mol_nums = []
            for entity_id in branch_list_entity_id:
                mol_nums.append(entity_number_of_molecules[int(entity_id) - 1]) # transforma o vetor de ids pra inteiro e retira 1 pq o id de entidades começa em 1 e o indice do vetor começa em 0

            olig_monosaccharide_dict = {"comp_id": branch_list_comp_id, "entry_id": mmcif_dict["_entry.id"][0], "oligossacaride": True, "entity_id": branch_list_entity_id, "comp_num": branch_list_comp_num, "mol_num": mol_nums}
            olig_monosaccharide_df = pd.DataFrame(data = olig_monosaccharide_dict)

        #non-oligosaccharide monosaccharides
        if "_pdbx_entity_nonpoly.entity_id" in mmcif_dict:

            nonpoly_entity_id = mmcif_dict["_pdbx_entity_nonpoly.entity_id"]
            nonpoly_entity_comp_id = mmcif_dict["_pdbx_entity_nonpoly.comp_id"]

            mol_nums = []
            for entity_id in nonpoly_entity_id:
                mol_nums.append(entity_number_of_molecules[int(entity_id) - 1]) # transforma o vetor de ids pra inteiro e retira 1 pq o id de entidades começa em 1 e o indice do vetor começa em 0

            monosaccharide_dict = {"comp_id": nonpoly_entity_comp_id, "entry_id": mmcif_dict["_entry.id"][0], "oligossacaride": False, "entity_id":  nonpoly_entity_id, "comp_num": None, "mol_num": mol_nums}
            monosaccharide_df = pd.DataFrame(data = monosaccharide_dict)

            carbo_dict = pd.read_csv("/home/douglas_lima/pdb/dicts/CCD_carbohydrate_list.tsv", sep = "\t", header = None, names = ['carbo_id', 'REF'])

            monosaccharide_df = monosaccharide_df[monosaccharide_df.comp_id.isin(carbo_dict["carbo_id"].values)]
           
        #junta os dataframes dos monossacarídeos
        monosaccharides = pd.concat([monosaccharides, olig_monosaccharide_df, monosaccharide_df], ignore_index=True)
        
    #mudar pra criar csv
    print(oligosaccharide_df)
    print(monosaccharides)


#fileNames = os.listdir("/home/douglas_lima/pdb/testesCif")
fileNames = ["/home/douglas_lima/pdb/testesCif/1v6u.cif", "/home/douglas_lima/pdb/testesCif/2wmg.cif", "/home/douglas_lima/pdb/testesCif/4of3.cif", "/home/douglas_lima/pdb/testesCif/5ebw.cif"]
separate(fileNames)

#print(monosaccharide_df_2)
#print(mmcif_dict["_entry.id"])
#/////////////////////////////////////////////////////////////////////////////
#
# os.chdir("/home/douglas_lima/pdb/dicts")
# comp_dict_file = open("/home/douglas_lima/pdb/dicts/components.cif")

# print(comp_dict_file)


# comp_dict = MMCIF2Dict(comp_dict_file.name)
# comp_dict_file.close()
#
#/////////////////////////////////////////////////////////////////////////////

#carbo_dict = pd.read_csv("/home/douglas_lima/pdb/dicts/CCD_carbohydrate_list.tsv", sep = "\t", header = None, names = ['carbo_id', 'REF'])

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