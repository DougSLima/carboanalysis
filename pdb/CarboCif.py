import os
import time
import re
from itertools import product
from itertools import repeat
import concurrent.futures
import pandas as pd
import warnings
from multiprocessing import Process, Pool
from pathlib import *
from Bio.PDB import *
from Bio.PDB.MMCIF2Dict import *
from Bio.PDB.PDBExceptions import PDBConstructionWarning 

warnings.simplefilter('ignore', PDBConstructionWarning) #ignorar warning (PDBConstructionWarning: WARNING: Chain B is discontinuous at line numeroDaLinha.)

# Remove da lista recebida todos os nomes de arquivos cujas estruturas apresentem resoluções maiores que a resolução máxima selecionada
def filter_maxResolution(fileNames, maxResolution):

    filteredFileNames = [] # lista de filtrados
    pdbParser = MMCIFParser() # parser Biopython

    for fileName in fileNames:
        file = open(fileName, 'r')
        structure = pdbParser.get_structure(file.name, file)
        file.close()
        
        if structure.header.get("resolution") <= maxResolution: # Testa se a resolução da estrutura é menor ou igual a resolução desejada (atributo da função: maxResolution)
            filteredFileNames.append(fileName)
    
    return filteredFileNames

# Remove da lista recebida todos os nomes de arquivos cujas estruturas apresentem OWAB maiores que o OWAB máximo selecionado
def filter_maxOWAB(fileNames, maxOWAB):

    filteredFileNames = [] # lista de filtrados
    pdbParser = FastMMCIFParser() # parser Biopython

    for fileName in fileNames:
        file = open(fileName, 'r')
        structure = pdbParser.get_structure(file.name, file)
        file.close()
        atomList = Selection.unfold_entities(structure, 'A')
        occupancy_x_bfactor = 0

        for atom in atomList:
            occupancy_x_bfactor += (atom.get_occupancy() * atom.get_bfactor())
        
        owab = occupancy_x_bfactor/len(atomList) # OWAB = média de ocupância * b_factor dos átomos da estrutura

        if owab <= maxOWAB:
            filteredFileNames.append(fileName)

    return filteredFileNames

# Remove da lista recebida todos os nomes de arquivos cujas estruturas apresentem métodos de estrutura diferentes do método escolhido
def filter_structureMethod(fileNames, structure_method):
    
    # Allowed values for structure_method:
    #
    # 'X-RAY DIFFRACTION'
    # 'ELECTRON CRYSTALLOGRAPHY'	
    # 'ELECTRON MICROSCOPY'	
    # 'EPR'	
    # 'FIBER DIFFRACTION'	
    # 'FLUORESCENCE TRANSFER'	
    # 'INFRARED SPECTROSCOPY'	
    # 'NEUTRON DIFFRACTION'	
    # 'POWDER DIFFRACTION'	
    # 'SOLID-STATE NMR'	
    # 'SOLUTION NMR'	
    # 'SOLUTION SCATTERING'
    # 'THEORETICAL MODEL'

    filteredFileNames = [] # lista de filtrados
    pdbParser = MMCIFParser()

    for fileName in fileNames:
        file = open(fileName, 'r')
        structure = pdbParser.get_structure(file.name, file)
        file.close()
        
        if structure.header.get("structure_method") == structure_method:
            filteredFileNames.append(fileName)
    
    return filteredFileNames

# Remove da lista recebida todos os nomes de arquivos cujas estruturas não apresentam carboidratos
def filter_containCarbo(fileNames):
    
    filteredFileNames = []
    filteredEntries = []
    carbo_dict = pd.read_csv("/home/douglas_lima/pdb/dicts/CCD_carbohydrate_list.tsv", sep = "\t", header = None, names = ['carbo_id', 'release_status']) # release status values: https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp.pdbx_release_status.html
    #carbo_dict = pd.read_csv("/home/douglas/carboanalysis/carboanalysis/pdb/dicts/CCD_carbohydrate_list.tsv", sep = "\t", header = None, names = ['carbo_id', 'release_status']) # release status values: https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp.pdbx_release_status.html
    for fileName in fileNames:
        
        mmcif_dict = MMCIF2Dict(fileName)
        chem_components = mmcif_dict["_chem_comp.id"]

        if(set(carbo_dict["carbo_id"].values) & set(chem_components)):
            
            filteredFileNames.append(fileName)
            filteredEntries.append(mmcif_dict["_entry.id"][0])

    # escreve um arquivo .txt com o nome das entradas cujas estrutuaras possuam carbohidratos
    with open("/home/douglas_lima/pdb/dataframes/carbo_entrys.txt", "w") as file:
    #with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/carbo_entrys.txt", "w") as file:
        for entry in filteredEntries:
            # escreve cada entry_id em uma nova linha
            file.write("%s\n" % entry)
        print('Done: carbo_entrys.txt')    

    return filteredFileNames

# Separa os monossacarideos e oligossacarideos em dois dataframes gerando dois arquivos .csv
# Chama: filter_containCarbo()
def separate(fileNames):
    
    fileNames = filter_containCarbo(fileNames)

    oligosaccharide_df = pd.DataFrame() # dataframe de oligossacarídeos
    olig_monosaccharide_df = pd.DataFrame() # dataframe de monossacarídeos contidos em oligossacarídeos
    monosaccharide_df = pd.DataFrame() # dataframe de monossacarídeos "isolados"
    olig_and_non_olig_monosaccharides = pd.DataFrame() # dataframe que junta os dois acima (monossacarídeos que pertencem a oligossacarídeos + monossacarídeos "isolados")
    monosaccharides = pd.DataFrame() # dataframe resultante após a adição de commom names e iupac names no dataframe olig_and_non_olig_monosaccharides

    for fileName in fileNames:

        # cria um dicionário a partir do arquivo .cif
        mmcif_dict = MMCIF2Dict(fileName)
        
        # entity data
        entity_ids = mmcif_dict["_entity.id"]
        entity_types = mmcif_dict["_entity.type"]
        entity_descriptions = mmcif_dict["_entity.pdbx_description"]
        entity_number_of_molecules = [eval(i) for i in mmcif_dict["_entity.pdbx_number_of_molecules"]] # usando eval() pra converter o numero de moleculas de str pra int
        entity_formula_weight = mmcif_dict["_entity.formula_weight"]

        # oligosaccharide DataFrame
        if "_pdbx_entity_branch.entity_id" in mmcif_dict:

            entity_dict = {"id": entity_ids, "entry_id": mmcif_dict["_entry.id"][0], "type": entity_types, "description": entity_descriptions,"mol_num": entity_number_of_molecules, "formula_weight": entity_formula_weight}
            entity_df = pd.DataFrame(data = entity_dict)
            
            branch_entity_id = mmcif_dict["_pdbx_entity_branch.entity_id"]
            branch_entity_type = mmcif_dict["_pdbx_entity_branch.type"]

            # usa a lista de ids de entidade e de tipo pra verificar quais entidades são oligossacarídeos
            olig_dic = dict(zip(branch_entity_id, branch_entity_type))

            entity_df = entity_df[entity_df.id.isin([k for k, v in olig_dic.items() if v == 'oligosaccharide'])]
            
            oligosaccharide_df = pd.concat([oligosaccharide_df, entity_df], ignore_index=True)

        # monosaccharide from oligosaccharide DataFrame
        if "_pdbx_entity_branch_list.entity_id" in mmcif_dict:

            branch_list_entity_id = mmcif_dict["_pdbx_entity_branch_list.entity_id"]
            branch_list_comp_id = mmcif_dict["_pdbx_entity_branch_list.comp_id"]
            branch_list_comp_num = mmcif_dict["_pdbx_entity_branch_list.num"]

            mol_nums = []
            for entity_id in branch_list_entity_id:
                mol_nums.append(entity_number_of_molecules[int(entity_id) - 1]) # transforma o vetor de ids pra inteiro e retira 1 pq o id de entidades começa em 1 e o indice do vetor começa em 0

            olig_monosaccharide_dict = {"comp_id": branch_list_comp_id, "entry_id": mmcif_dict["_entry.id"][0], "oligossacaride": True, "entity_id": branch_list_entity_id, "comp_num": branch_list_comp_num, "mol_num": mol_nums}
            olig_monosaccharide_df = pd.DataFrame(data = olig_monosaccharide_dict)

        # non-oligosaccharide monosaccharides DataFrame
        if "_pdbx_entity_nonpoly.entity_id" in mmcif_dict:

            nonpoly_entity_id = mmcif_dict["_pdbx_entity_nonpoly.entity_id"]
            nonpoly_entity_comp_id = mmcif_dict["_pdbx_entity_nonpoly.comp_id"]

            mol_nums = []
            for entity_id in nonpoly_entity_id:
                mol_nums.append(entity_number_of_molecules[int(entity_id) - 1]) # transforma o vetor de ids pra inteiro e retira 1 pq o id de entidades começa em 1 e o indice do vetor começa em 0

            monosaccharide_dict = {"comp_id": nonpoly_entity_comp_id, "entry_id": mmcif_dict["_entry.id"][0], "oligossacaride": False, "entity_id":  nonpoly_entity_id, "comp_num": None, "mol_num": mol_nums}
            monosaccharide_df = pd.DataFrame(data = monosaccharide_dict)

            #carbo_dict = pd.read_csv("/home/douglas_lima/pdb/dicts/CCD_carbohydrate_list.tsv", sep = "\t", header = None, names = ['carbo_id', 'REF'])
            carbo_dict = pd.read_csv("/home/douglas/carboanalysis/carboanalysis/pdb/dicts/CCD_carbohydrate_list.tsv", sep = "\t", header = None, names = ['carbo_id', 'REF'])

            monosaccharide_df = monosaccharide_df[monosaccharide_df.comp_id.isin(carbo_dict["carbo_id"].values)]
           
        # junta os DataFrames dos monossacarídeos que pertencem com os que não pertencem a oligossacarídeos
        olig_and_non_olig_monosaccharides = pd.concat([olig_monosaccharide_df, monosaccharide_df], ignore_index=True) 
        
    # chem comp identifier   
        identifier_comp_id = mmcif_dict["_pdbx_chem_comp_identifier.comp_id"]
        identifier_type = mmcif_dict["_pdbx_chem_comp_identifier.type"]
        identifier_identifier = mmcif_dict["_pdbx_chem_comp_identifier.identifier"]

        identifier_dict = {"comp_id": identifier_comp_id, "type": identifier_type, "identifier": identifier_identifier}
        identifier_df = pd.DataFrame(data=identifier_dict)

        commom_names = []
        iupac_symbols = []

        for comp_id in olig_and_non_olig_monosaccharides.comp_id:
            commom_names.append(identifier_df[(identifier_df.comp_id == comp_id) & (identifier_df.type == 'COMMON NAME')]["identifier"].values[0]) 
            iupac_symbols.append(identifier_df[(identifier_df.comp_id == comp_id) & (identifier_df.type == 'IUPAC CARBOHYDRATE SYMBOL')]["identifier"].values[0]) 
    
        olig_and_non_olig_monosaccharides["commom_name"] = commom_names
        olig_and_non_olig_monosaccharides["iupac_symbol"] = iupac_symbols      
        
        monosaccharides = pd.concat([monosaccharides, olig_and_non_olig_monosaccharides], ignore_index=True)
    
    #monosaccharides.to_csv(path_or_buf="/home/douglas_lima/pdb/dataframes/monosaccharides.csv") # escreve o .csv a partir do DataFrame dos monossacarídeos
    #oligosaccharide_df.to_csv(path_or_buf="/home/douglas_lima/pdb/dataframes/oligosaccharides.csv") # escreve o .csv a partir do DataFrame dos oligossacarídeos
    
    monosaccharides.to_csv(path_or_buf="/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/monosaccharides.csv") # escreve o .csv a partir do DataFrame dos monossacarídeos
    oligosaccharide_df.to_csv(path_or_buf="/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/oligosaccharides.csv") # escreve o .csv a partir do DataFrame dos oligossacarídeos

    # monosaccharides_sum DataFrame (.csv)
    comp_ids = []
    sums = []
    commom_names = []
    iupac_symbols = []

    # monta as listas das variáveis pra formar o DataFrame
    for carbo in monosaccharides.comp_id.unique(): # realiza a iteração para cada TAG (comp_id) única 
        comp_ids.append(carbo)
        sums.append(monosaccharides.loc[monosaccharides['comp_id'] == carbo, 'mol_num'].sum()) # soma a coluna do número de moléculas pra todas as linhas que tiverem a mesma TAG de monômeros
        commom_names.append(monosaccharides.loc[monosaccharides['comp_id'] == carbo, 'commom_name'].unique()[0])
        iupac_symbols.append(monosaccharides.loc[monosaccharides['comp_id'] == carbo, 'iupac_symbol'].unique()[0])

    carbo_dict = {"comp_id": comp_ids, "sum": sums, "commom_name": commom_names, "iupac_symbol": iupac_symbols}
    carbo_df = pd.DataFrame(data = carbo_dict)
    #carbo_df.to_csv(path_or_buf="/home/douglas_lima/pdb/dataframes/monosaccharides_sum.csv") # escreve o .csv com a soma da ocorrência dos monossacarídeos
    carbo_df.to_csv(path_or_buf="/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/monosaccharides_sum.csv")

def bfactor_values(fileName):

    mmcif_dict = MMCIF2Dict(fileName)#_atom_site.B_iso_or_equiv
    
    atom_site_id = mmcif_dict['_atom_site.id'] 
    atom_group = mmcif_dict['_atom_site.group_PDB']
    entity_id = mmcif_dict['_atom_site.label_entity_id']
    bfactors = mmcif_dict['_atom_site.B_iso_or_equiv']
    
    

def bfactor_parser(fileNames):

    pdbParser = MMCIFParser() # parser Biopython

    for fileName in fileNames:
        file = open(fileName, 'r')
        structure = pdbParser.get_structure(file.name, file)
        file.close()
        atomList = Selection.unfold_entities(structure, 'A')
        a = atomList[0]
        print(a.get_bfactor())
        a = atomList[35]
        print(a.get_bfactor())
        a = atomList[1200]
        print(a.get_bfactor())

os.chdir("/home/douglas_lima/pdb/testesCif")
fileNames = os.listdir("/home/douglas_lima/pdb/testesCif")
# os.chdir("/home/douglas/carboanalysis/data/unzipped")
# fileNames = os.listdir("/home/douglas/carboanalysis/data/unzipped")

# filtrados = filter_structureMethod(fileNames, 'X-RAY DIFFRACTION')
# filtrados = filter_containCarbo(filtrados)


# print(fileNames)
# print(filtrados)

filtrados = ["2wmg.cif"]
bfactor_parser(filtrados)
with Pool() as pool:
        print("bfactor: Starting...")
        results = [i for i in pool.map(bfactor_values, filtrados) if i is not None]
        print("bfactor: Done!")