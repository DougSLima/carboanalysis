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

# Remove da lista recebida todos os nomes de arquivos cujas estruturas apresentem métodos de estrutura diferentes do método escolhido, escrevendo em um arquivo de saída
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
            # escreve um arquivo .txt com o nome das entradas
            with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/xray_diffraction_carbo_entrys.txt", "a") as file:
                file.write("%s\n" % fileName) 
            return fileName

#Escreve as resoluções dos arquivos
def write_resolution(fileName):

    mmcif_dict = MMCIF2Dict(fileName)
    
    try:
        resolution = float(mmcif_dict["_refine.ls_d_res_high"][0])
    except ValueError as error:
        with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/SIC/valueerror_entrys.txt", "a") as file:
            file.write("%s\n" % fileName)
        return None
    except KeyError:
        with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/SIC/keyerror_entrys.txt", "a") as file:
            file.write("%s\n" % fileName) 
        return None
    
    try:
        if("_refine.ls_d_res_high" in mmcif_dict.keys()):
            res_dict = {"entry": mmcif_dict['_entry.id'], "resolution": resolution}
            res_df = pd.DataFrame(data = res_dict)
            res_df.to_csv(path_or_buf="/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/SIC/res.csv", mode='a', index=False, header=False, sep=";")
    except KeyError:
        with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/SIC/keyerror_entrys.txt", "a") as file:
            file.write("%s\n" % fileName) 
        return None

#Separa os monômeros
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
        # !!!!!!!!!!!!!!!!! PROBLEMA AQUI !!!!!!!!!!!!!!!!!!
        #oligossaccharide_df.to_csv(path_or_buf="/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/SIC/oligossaccharides.csv", mode='a', index=False, header=False) # escreve o .csv de oligossacarídeos

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
    
    olig_and_non_olig_monossaccharides.to_csv(path_or_buf="/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/SIC/monossaccharides.csv", mode='a', index=False, header=False) # escreve o .csv de monossacarídeos
    # oligossaccharide_df.to_csv(path_or_buf="/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/oligossaccharides.csv") # escreve o .csv de oligossacarídeos

def filter_maxResolution(fileName):

    mmcif_dict = MMCIF2Dict(fileName)
    
    try:
        resolution = float(mmcif_dict["_refine.ls_d_res_high"][0])
    except ValueError as error:
        #with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/valueerror_entrys.txt", "a") as file:
            #file.write("%s\n" % fileName)
        return None
    except KeyError:
        #with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/keyerror_entrys.txt", "a") as file:
            #file.write("%s\n" % fileName) 
        return None
    
    try:
        if("_refine.ls_d_res_high" in mmcif_dict.keys()):
            if(resolution <= 2.0):
                return fileName
    except KeyError:
        #with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/keyerror_entrys.txt", "a") as file:
            #file.write("%s\n" % fileName) 
        return None
    
def filter_OWAB(fileName):
    
    print("OWAB filtering: " + fileName)

    mmcif_dict = MMCIF2Dict(fileName)

    #pandas df
    atom_df_dict = {"bfactor": mmcif_dict['_atom_site.B_iso_or_equiv'], "occupancy": mmcif_dict['_atom_site.occupancy']}
    atom_df = pd.DataFrame(data = atom_df_dict)

    #string to float
    atom_df['bfactor'] = atom_df["bfactor"].astype(float)
    atom_df['occupancy'] = atom_df["occupancy"].astype(float)
    
    atom_df['m'] = atom_df['bfactor'] * atom_df['occupancy']

    if(atom_df['m'].mean() <= 60):
        return fileName
    else:
        return None
    
def write_carbo(fileName):
    print("Writing: " + fileName)
    # escreve um arquivo .txt com o nome das entradas cujas estrutuaras possuam carbohidratos
    with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/SIC/SIC2023_carbo_entrys_res_owab_filtered.txt", "a") as file:
        file.write("%s\n" % fileName) 

#Calcula a diferença do bfactor entre carboidrato e proteina (salva somente a maior)
def bfactor_values(fileName):

    print("B - factor values: " + fileName)

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

    #print(hetatm_df)
    #remove moléculas de água
    hetatm_df = hetatm_df.loc[hetatm_df['comp'] != 'HOH']

    #print(hetatm_df)

    #Separa so os carboidratos
    carbo_dict = pd.read_csv("/home/douglas/carboanalysis/carboanalysis/pdb/dicts/CCD_carbohydrate_list.tsv", sep = "\t", header = None, names = ['carbo_id', 'release_status'])
    # carbo_dict = pd.read_csv("/home/douglas_lima/pdb/dicts/CCD_carbohydrate_list.tsv", sep = "\t", header = None, names = ['carbo_id', 'release_status'])
    carbo_list = carbo_dict["carbo_id"].values

    hetatm_df = hetatm_df.loc[hetatm_df['comp'].isin(carbo_list)]

    #print(hetatm_df)

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
            
                # print(comp + ' trocou para: ' + row['comp'])
                comp = row["comp"]
                # print(bfactors)
                # print("Media: " + str(media))
                bfactors = []
                # print("Maior media: " + str(maiorMedia))
                # print("Maior comp: " + maiorComp)
            
            bfactors.append(row["bfactors"])

    entry_dict = {"entry": mmcif_dict['_entry.id'], "polymer_mean": polymer_mean, "mbfctor_comp": maiorComp, "mbfcator_mean": maiorMedia, "diff": polymer_mean - maiorMedia}
    entry_df = pd.DataFrame(data = entry_dict)
    entry_df.to_csv(path_or_buf="/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/SIC/bfactors.csv", mode='a', index=False, header=False, sep=";")

if __name__ == '__main__':

    start = time.time()
    #certin
    #df = pd.read_csv("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/xray_diffraction_carbo_entrys.txt", names = ['entry_filename'])

    #df = pd.read_csv("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/carbo_entrys_res_owab_filtered_nd.txt", names = ['entry_filename'])
    df = pd.read_csv("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/SIC/SIC2023_carbo_entrys_res_owab_filtered.txt", names = ['entry_filename'])

    fileNames = df['entry_filename'].values
    
    os.chdir("/home/douglas/carboanalysis/data/unzipped")

    with Pool() as pool:
        # print("Writing xray diff...")
        # results = [i for i in pool.map(filter_structureMethod, fileNames) if i is not None]
        # print("Done!")
        # print("Writing entries resolutions...")
        # results = [i for i in pool.map(write_resolution, fileNames) if i is not None]
        # print("Writing entries resolutions: Done!")
        # print("Sep...")
        # results = [i for i in pool.map(separate, fileNames) if i is not None]
        # print("Sep: Done!")
        # print("Res Filter: ...")
        # results = [i for i in pool.map(filter_maxResolution, fileNames) if i is not None]
        # print("Res filter: Done!")
        # print("Owab Filter: ...")
        # results = [i for i in pool.map(filter_OWAB, results) if i is not None]
        # print("Owab Filter: Done!")
        # print("Writing carbo_entrys file...")
        # results = [i for i in pool.map(write_carbo, results) if i is not None]
        # print("Writing Done!")
        print("Bfactoring...")
        results = [i for i in pool.map(bfactor_values, fileNames) if i is not None]
        print("Bfactor parsing Done!")


    print("thread time: ", time.time() - start)