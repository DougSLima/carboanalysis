import os
import time
import re
from itertools import product
from itertools import repeat
import concurrent.futures
import pandas as pd
import statistics
import warnings
from multiprocessing import Process, Pool
from pathlib import *
from Bio.PDB import *
from Bio.PDB.MMCIF2Dict import *
from Bio.PDB.PDBExceptions import PDBConstructionWarning 

warnings.simplefilter('ignore', PDBConstructionWarning) #ignorar warning (PDBConstructionWarning: WARNING: Chain B is discontinuous at line numeroDaLinha.)

# Remove da lista recebida todos os nomes de arquivos cujas estruturas apresentem resoluções maiores que a resolução máxima selecionada
def filter_maxResolution(fileName):

    mmcif_dict = MMCIF2Dict(fileName)

    if "_reflns.d_resolution_high" in mmcif_dict.keys():
        if float(mmcif_dict["_reflns.d_resolution_high"][0]) <= 2.5:
            return fileName

# Remove da lista recebida todos os nomes de arquivos cujas estruturas apresentem métodos de estrutura diferentes do método escolhido
def filter_structureMethod(fileName):
    
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

    mmcif_dict = MMCIF2Dict(fileName)

    if '_exptl.method' in mmcif_dict.keys():
        if mmcif_dict['_exptl.method'][0] == 'X-RAY DIFFRACTION':
            return fileName

# Remove da lista recebida todos os nomes de arquivos cujas estruturas não apresentam carboidratos
def filter_containCarbo(fileName):
    
    print("Analysing: " + fileName)
    carbo_dict = pd.read_csv("/home/douglas/carboanalysis/carboanalysis/pdb/dicts/CCD_carbohydrate_list.tsv", sep = "\t", header = None, names = ['carbo_id', 'release_status']) # release status values: https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp.pdbx_release_status.html
    #carbo_dict = pd.read_csv("/home/douglas/carboanalysis/carboanalysis/pdb/dicts/CCD_carbohydrate_list.tsv", sep = "\t", header = None, names = ['carbo_id', 'release_status']) # release status values: https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp.pdbx_release_status.html
    
    mmcif_dict = MMCIF2Dict(fileName)
    chem_components = mmcif_dict["_chem_comp.id"]

    if(set(carbo_dict["carbo_id"].values) & set(chem_components)):
        return fileName

def write_carbo(fileName):
    print("Writing: " + fileName)
    # escreve um arquivo .txt com o nome das entradas cujas estrutuaras possuam carbohidratos
    with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/carbo_entrys.txt", "a") as file:
    #with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/carbo_entrys.txt", "a") as file:
        file.write("%s\n" % fileName) 

# Separa os monossacarideos e oligossacarideos em dois dataframes gerando dois arquivos .csv
def separate(fileName):

    oligossaccharide_df = pd.DataFrame() # dataframe de oligossacarídeos
    olig_monosaccharide_df = pd.DataFrame() # dataframe de monossacarídeos contidos em oligossacarídeos
    monosaccharide_df = pd.DataFrame() # dataframe de monossacarídeos "isolados"
    olig_and_non_olig_monossaccharides = pd.DataFrame() # dataframe que junta os dois acima (monossacarídeos que pertencem a oligossacarídeos + monossacarídeos "isolados")
    monossaccharides = pd.DataFrame() # dataframe resultante após a adição de commom names e iupac names no dataframe olig_and_non_olig_monossaccharides

    print("Separating: " + fileName)

    entity_df = pd.DataFrame()
    olig_monosaccharide_df = pd.DataFrame() # dataframe de monossacarídeos contidos em oligossacarídeos
    monosaccharide_df = pd.DataFrame() # dataframe de monossacarídeos "isolados"

    # cria um dicionário a partir do arquivo .cif
    mmcif_dict = MMCIF2Dict(fileName)
    
    # entity data
    entity_ids = mmcif_dict["_entity.id"]
    entity_types = mmcif_dict["_entity.type"]
    entity_descriptions = mmcif_dict["_entity.pdbx_description"]
    entity_number_of_molecules = [eval(i) for i in mmcif_dict["_entity.pdbx_number_of_molecules"]] # usando eval() pra converter o numero de moleculas de str pra int
    entity_formula_weight = mmcif_dict["_entity.formula_weight"]
    
    # oligossaccharide DataFrame
    if "_pdbx_entity_branch.entity_id" in mmcif_dict.keys():

        entity_dict = {"id": entity_ids, "entry_id": mmcif_dict["_entry.id"][0], "type": entity_types, "description": entity_descriptions,"mol_num": entity_number_of_molecules, "formula_weight": entity_formula_weight}
        entity_df = pd.DataFrame(data = entity_dict)
        
        branch_entity_id = mmcif_dict["_pdbx_entity_branch.entity_id"]
        branch_entity_type = mmcif_dict["_pdbx_entity_branch.type"]

        # usa a lista de ids de entidade e de tipo pra verificar quais entidades são oligossacarídeos
        olig_dic = dict(zip(branch_entity_id, branch_entity_type))

        oligossaccharide_df = entity_df[entity_df.id.isin([k for k, v in olig_dic.items() if v == 'oligossaccharide'])]
        
        oligossaccharide_df.to_csv(path_or_buf="/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/oligossaccharides.csv", mode='a', index=False, header=False) # escreve o .csv de oligossacarídeos

    # monosaccharide from oligossaccharide DataFrame
    if "_pdbx_entity_branch_list.entity_id" in mmcif_dict.keys():

        branch_list_entity_id = mmcif_dict["_pdbx_entity_branch_list.entity_id"]
        branch_list_comp_id = mmcif_dict["_pdbx_entity_branch_list.comp_id"]
        branch_list_comp_num = mmcif_dict["_pdbx_entity_branch_list.num"]

        mol_nums = []
        for entity_id in branch_list_entity_id:
            mol_nums.append(entity_number_of_molecules[int(entity_id) - 1]) # transforma o vetor de ids pra inteiro e retira 1 pq o id de entidades começa em 1 e o indice do vetor começa em 0
        
        olig_monosaccharide_dict = {"comp_id": branch_list_comp_id, "entry_id": mmcif_dict["_entry.id"][0], "oligossaccharide": True, "entity_id": branch_list_entity_id, "comp_num": branch_list_comp_num, "mol_num": mol_nums}
        olig_monosaccharide_df = pd.DataFrame(data = olig_monosaccharide_dict)

    # non-oligossaccharide monossaccharides DataFrame
    if "_pdbx_entity_nonpoly.entity_id" in mmcif_dict.keys():

        nonpoly_entity_id = mmcif_dict["_pdbx_entity_nonpoly.entity_id"]
        nonpoly_entity_comp_id = mmcif_dict["_pdbx_entity_nonpoly.comp_id"]

        mol_nums = []
        for entity_id in nonpoly_entity_id:
            mol_nums.append(entity_number_of_molecules[int(entity_id) - 1]) # transforma o vetor de ids pra inteiro e retira 1 pq o id de entidades começa em 1 e o indice do vetor começa em 0

        monosaccharide_dict = {"comp_id": nonpoly_entity_comp_id, "entry_id": mmcif_dict["_entry.id"][0], "oligossaccharide": False, "entity_id":  nonpoly_entity_id, "comp_num": None, "mol_num": mol_nums}
        monosaccharide_df = pd.DataFrame(data = monosaccharide_dict)

        carbo_dict = pd.read_csv("/home/douglas/carboanalysis/carboanalysis/pdb/dicts/CCD_carbohydrate_list.tsv", sep = "\t", header = None, names = ['carbo_id', 'REF'])

        monosaccharide_df = monosaccharide_df[monosaccharide_df.comp_id.isin(carbo_dict["carbo_id"].values)]
        
    # junta os dataframes dos monossacarídeos
    olig_and_non_olig_monossaccharides = pd.concat([olig_monosaccharide_df, monosaccharide_df], ignore_index=True)

    chem_comp_dict = {"comp_id": mmcif_dict['_chem_comp.id'], "linking": mmcif_dict['_chem_comp.type'], "name": mmcif_dict['_chem_comp.name']}
    chem_comp_df = pd.DataFrame(data = chem_comp_dict)

    names = []
    
    for comp_id in olig_and_non_olig_monossaccharides['comp_id'].values:
        names.append(chem_comp_df[chem_comp_df.comp_id == comp_id]["name"].values[0])    
    
    olig_and_non_olig_monossaccharides["name"] = names
    
    # chem comp identifier
    # identifier_comp_id = mmcif_dict["_pdbx_chem_comp_identifier.comp_id"]
    # identifier_type = mmcif_dict["_pdbx_chem_comp_identifier.type"]
    # identifier_identifier = mmcif_dict["_pdbx_chem_comp_identifier.identifier"]
    # print(mmcif_dict["_pdbx_chem_comp_identifier.type"])
    # print(type(mmcif_dict["_pdbx_chem_comp_identifier.type"]))

    # identifier_dict = {"comp_id": identifier_comp_id, "type": identifier_type, "identifier": identifier_identifier}
    # identifier_df = pd.DataFrame(data=identifier_dict)

    # commom_names = []
    # iupac_symbols = []

    # for comp_id in olig_and_non_olig_monossaccharides.comp_id:
    #     commom_names.append(identifier_df[(identifier_df.comp_id == comp_id) & (identifier_df.type == 'COMMON NAME')]["identifier"].values) 
    #     iupac_symbols.append(identifier_df[(identifier_df.comp_id == comp_id) & (identifier_df.type == 'IUPAC CARBOHYDRATE SYMBOL')]["identifier"].values) 

    # olig_and_non_olig_monossaccharides["commom_name"] = commom_names
    # olig_and_non_olig_monossaccharides["iupac_symbol"] = iupac_symbols      
    
    #monossaccharides = pd.concat([monossaccharides, olig_and_non_olig_monossaccharides], ignore_index=True)
    
    olig_and_non_olig_monossaccharides.to_csv(path_or_buf="/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/monossaccharides.csv", mode='a', index=False, header=False) # escreve o .csv de monossacarídeos
    # oligossaccharide_df.to_csv(path_or_buf="/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/oligossaccharides.csv") # escreve o .csv de oligossacarídeos

def unique_monossaccharides():

    monossaccharides = pd.read_csv('/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/monossaccharides.csv', index_col=None, names=['comp_id', 'entry_id', 'oligossaccharide', 'entity_id', 'comp_num', 'mol_num', 'name'], header=None)

    # monossaccharides_sum DataFrame (.csv)
    comp_ids = []
    sums = []
    names = []

    for carbo in monossaccharides.comp_id.unique(): # realiza a iteração para cada TAG (comp_id) única 
        comp_ids.append(carbo)
        sums.append(monossaccharides.loc[monossaccharides['comp_id'] == carbo, 'mol_num'].sum())
        names.append(monossaccharides.loc[monossaccharides['comp_id'] == carbo, 'name'].unique()[0])

    carbo_dict = {"comp_id": comp_ids, "sum": sums, "name": names}
    carbo_df = pd.DataFrame(data = carbo_dict)
    carbo_df.sort_values(by=['sum'], ascending=False)
    carbo_df.to_csv(path_or_buf="/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/monossaccharides_sum.csv") # escreve o .csv com a soma da ocorrência dos monossacarídeos

def bfactor_values(fileName):

    mmcif_dict = MMCIF2Dict(fileName)#_atom_site.B_iso_or_equiv
    
    # atom_group = mmcif_dict['_atom_site.group_PDB']
    # atom_site_id = mmcif_dict['_atom_site.id'] 
    # atom_comp_id = mmcif_dict['_atom_site.label_comp_id']
    # atom_label_id = mmcif_dict['_atom_site.label_atom_id']
    # bfactors = mmcif_dict['_atom_site.B_iso_or_equiv']
    # entity_id = mmcif_dict['_atom_site.label_entity_id']
    # entity_seq_num = mmcif_dict['_atom_site.label_seq_id']
    
    # junta os dataframes dos monossacarídeos
    atom_df_dict = {"group": mmcif_dict['_atom_site.group_PDB'], "atom_id": mmcif_dict['_atom_site.id'], "comp":  mmcif_dict['_atom_site.label_comp_id'], "atom_symbol": mmcif_dict["_atom_site.type_symbol"], "atom_label": mmcif_dict['_atom_site.label_atom_id'], "bfactors": mmcif_dict['_atom_site.B_iso_or_equiv'], "entity_id": mmcif_dict['_atom_site.label_entity_id'],
                    "entity_seq_num": mmcif_dict['_atom_site.label_seq_id']}
    atom_df = pd.DataFrame(data = atom_df_dict)

    atom_df['bfactors'] = atom_df["bfactors"].astype(float)
    #pega oligossacarideos se existem e mantém somente os ids unicos

    polymer_mean = atom_df.loc[atom_df['group'] == 'ATOM']["bfactors"].mean()
    #polymer_median = atom_df.loc[atom_df['group'] == 'ATOM']["bfactors"].median()

    #separa os hetero atomos
    hetatm_df = atom_df.loc[atom_df['group'] == 'HETATM']

    print(hetatm_df)
    #remove moléculas de água
    hetatm_df = hetatm_df.loc[hetatm_df['comp'] != 'HOH']

    print(hetatm_df)

    #Separa so os carboidratos
    #carbo_dict = pd.read_csv("/home/douglas/carboanalysis/carboanalysis/pdb/dicts/CCD_carbohydrate_list.tsv", sep = "\t", header = None, names = ['carbo_id', 'release_status'])
    carbo_dict = pd.read_csv("/home/douglas_lima/pdb/dicts/CCD_carbohydrate_list.tsv", sep = "\t", header = None, names = ['carbo_id', 'release_status'])
    carbo_list = carbo_dict["carbo_id"].values

    hetatm_df = hetatm_df.loc[hetatm_df['comp'].isin(carbo_list)]

    print(hetatm_df)

    #reseta os index
    hetatm_df = hetatm_df.reset_index()

    bfactors = []
    maiorMedia = 0
    maiorComp = ""
    comp = ""

    for index, row in hetatm_df.iterrows():
        
        if(index == 0):
            comp = row["comp"]
            bfactors.append(row["bfactors"])          
        else:
            if(row["comp"] != comp or row["atom_label"] == 'C1' or index == len(hetatm_df.index) - 1):
                
                if(index == len(hetatm_df.index) - 1):
                    bfactors.append(row["bfactors"])

                media = statistics.mean(bfactors)
                #mediana = statistics.median(bfactors)

                if(media > maiorMedia):
                    maiorMedia = media
                    maiorComp = comp
            
                print(comp + ' trocou para: ' + row['comp'])
                comp = row["comp"]
                print(bfactors)
                print("Media: " + str(media))
                bfactors = []
                print("Maior media: " + str(maiorMedia))
                print("Maior comp: " + maiorComp)
            
            bfactors.append(row["bfactors"])
    
    print("Media polimero: ", polymer_mean)
    print("Maior media: ", maiorMedia)
    print("maior comp:", maiorComp)

    entry_dict = {"entry": mmcif_dict['_entry.id'], "polymer_mean": polymer_mean, "mbfctor_comp": maiorComp, "mbfcator_mean": maiorMedia, "diff": polymer_mean - maiorMedia}
    entry_df = pd.DataFrame(data = entry_dict)
    entry_df.to_csv(path_or_buf="/home/douglas_lima/pdb/dataframes/bfactors.csv", mode='a', index=False, header=False, sep=";")



def bfactor(fileName):
    p = MMCIFParser()
    structure = p.get_structure(fileName, fileName)
    residues = structure.get_residues()

    for r in residues:
        print(r)

    






os.chdir("/home/douglas_lima/pdb/testesBfactor")
fileNames = os.listdir("/home/douglas_lima/pdb/testesBfactor")
# os.chdir("/home/douglas/carboanalysis/data/unzipped")
# fileNames = os.listdir("/home/douglas/carboanalysis/data/unzipped")

# filtrados = filter_structureMethod(fileNames, 'X-RAY DIFFRACTION')
# filtrados = filter_containCarbo(filtrados)


# print(fileNames)
# print(filtrados)

filtrados = ["2wmg.cif"]
# for i in filtrados:
#     bfactor_values(i)
with Pool() as pool:
        print("bfactor: Starting...")
        results = [i for i in pool.map(bfactor_values, fileNames) if i is not None]
        print("bfactor: Done!")

# if __name__ == '__main__':

#     start = time.time()
    
#     print("Unique processing...")

#     unique_monossaccharides()

#     print("Done!")
   
#     print("thread time: ", time.time() - start)
    