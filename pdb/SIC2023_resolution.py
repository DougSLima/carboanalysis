import os
import subprocess
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
pd.options.mode.chained_assignment = None 

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
    atom_df_dict = {"group": mmcif_dict['_atom_site.group_PDB'], "atom_id": mmcif_dict['_atom_site.id'], 
                    "comp":  mmcif_dict['_atom_site.label_comp_id'], "atom_symbol": mmcif_dict["_atom_site.type_symbol"], 
                    "atom_label": mmcif_dict['_atom_site.label_atom_id'], "bfactors": mmcif_dict['_atom_site.B_iso_or_equiv'], "entity_id": mmcif_dict['_atom_site.label_entity_id'],
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

#Coleta os dados das ligações glicosídicas
def find_linkages(fileName):
    
    #Cria um dicionário a partir do arquivo .cif
    mmcif_dict = MMCIF2Dict(fileName)

    try:
        print('Linking: ' + fileName)

        #Coleta informações das Branched entities
        branch_dict = {"entity_id": mmcif_dict['_entity.id'], 
                        "num_of_molecules": mmcif_dict['_entity.pdbx_number_of_molecules']}

        #Transforma num dataframe pandas
        branch_df = pd.DataFrame(data = branch_dict)

        #converte a coluna "num_of_molecules" para inteiro
        branch_df['num_of_molecules'] = branch_df['num_of_molecules'].astype(int)

        #Coleta informações das ligações
        linkage_dict = {"entry_id": mmcif_dict['_entry.id'],
                        "link_id": mmcif_dict['_pdbx_entity_branch_link.link_id'], 
                        "entity_id": mmcif_dict['_pdbx_entity_branch_link.entity_id'], 
                        "branch_1_id":  mmcif_dict['_pdbx_entity_branch_link.entity_branch_list_num_1'], 
                        "comp_1_id": mmcif_dict["_pdbx_entity_branch_link.comp_id_1"], 
                        "atom_1_id": mmcif_dict['_pdbx_entity_branch_link.atom_id_1'], 
                        "leaving_atom_1_id": mmcif_dict['_pdbx_entity_branch_link.leaving_atom_id_1'],
                        "branch_2_id":  mmcif_dict['_pdbx_entity_branch_link.entity_branch_list_num_2'], 
                        "comp_2_id": mmcif_dict["_pdbx_entity_branch_link.comp_id_2"], 
                        "atom_2_id": mmcif_dict['_pdbx_entity_branch_link.atom_id_2'], 
                        "leaving_atom_2_id": mmcif_dict['_pdbx_entity_branch_link.leaving_atom_id_2'],
                        "order": mmcif_dict['_pdbx_entity_branch_link.value_order']}
        
        #Transforma num dataframe pandas
        linkage_df = pd.DataFrame(data = linkage_dict)
        
        #Coleta informações de nomenclatura dos açúcares
        identifier_dict = {"comp_id": mmcif_dict['_pdbx_chem_comp_identifier.comp_id'], 
                        "identifier_type": mmcif_dict['_pdbx_chem_comp_identifier.type'], 
                        "identifier":  mmcif_dict['_pdbx_chem_comp_identifier.identifier']}
        
        #Transforma num dataframe pandas
        identifier_df = pd.DataFrame(data = identifier_dict)

        #Leva em consideração o número de moleculas de cada entidade
        num_of_molecules_list = []
        for index, row in linkage_df.iterrows():
            num_of_molecules_list.append(branch_df.loc[branch_df['entity_id'] == row['entity_id'], 'num_of_molecules'].values[0])
        
        #Adiciona esse número como uma nova coluna do dataframe linkage_df
        linkage_df['num_of_molecules'] = num_of_molecules_list
        

        #Escreve a informação das ligações num arquivo .csv
        linkage_df.to_csv(path_or_buf="/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/SIC/all_linkages_v2.csv", mode='a', index=False, header=False, sep=";") 

    except ValueError as error:
        return None
    except KeyError as error:
        return None

def write_sugar_line_pdb(fileName, row, sugar, sugar_first_id):
    try:
        with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/pdb_sugars/" + fileName[:4] + "_" + sugar + "_" + sugar_first_id + ".pdb", "a") as f:
            
            # Formatar a linha de acordo com o padrão PDB
            pdb_line = "{:<6}{:>5} {:<4} {:<3} {:>1}{:>4}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:<2}".format(
                row['group'], row['id'], row['label_atom_id'], row['label_comp_id'], row['label_asym_id'], row['label_entity_id'],
                row['Cartn_x'], row['Cartn_y'], row['Cartn_z'], row['occupancy'], row['B_iso_or_equiv'],
                row['type_symbol']
            )
            f.write(pdb_line + '\n')
            
    except Exception as e:
        with open('/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/logs/log_de_erros.txt', 'a') as f:
            f.write(f"Erro: {str(e)}\n" + " file: " + fileName)
        print(row)

def is_piranose_or_furanose(sugar_id, mmcif_dict):
    chem_comp_dict = {"comp_id": mmcif_dict['_chem_comp.id'], "name": mmcif_dict['_chem_comp.name']}
    chem_comp_df = pd.DataFrame(data = chem_comp_dict)

    if("pyranose" in chem_comp_df[chem_comp_df.comp_id == sugar_id]["name"].values[0]):
        print("pyr")
        return "pyranose"

    elif("furanose" in chem_comp_df[chem_comp_df.comp_id == sugar_id]["name"].values[0]):
        print("fur")
        return "furanose"
    
def alter_dat(atoms_dict, ring_type, fileName, sugar, sugar_first_id):
    print("entro no alter")
    colvar_name = "/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/pcukering/colvar/" + fileName[:4] + "_" + sugar + "_" + sugar_first_id

    if(ring_type == "pyranose"):
        print("foi pir")
        try:
            o5 = atoms_dict['O5']
            c1 = atoms_dict['C1']
            c2 = atoms_dict['C2']
            c3 = atoms_dict['C3']
            c4 = atoms_dict['C4']
            c5 = atoms_dict['C5']
            print(atoms_dict)
            with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/puckering/puck.dat", "w") as f:
                f.write("#plumed.dat\n" +
                        "puck: PUCKERING ATOMS=" + 
                        o5 + ',' +
                        c1 + ',' +
                        c2 + ',' +
                        c3 + ',' +
                        c4 + ',' +
                        c5 + "\n" + 
                        "PRINT ARG=puck.* FILE=" + colvar_name)

        except Exception as e:
            with open('/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/logs/log_de_erros.txt', 'a') as f:
                f.write(f"Erro: {str(e)}\n" + " file: " + ring_type + " - " + atoms_dict)

    elif(ring_type == "furanose"):
        try:
            c4 = atoms_dict['C4']
            o4 = atoms_dict['O4']
            c1 = atoms_dict['C1']
            c2 = atoms_dict['C2']
            c3 = atoms_dict['C3']

            with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/puckering/puck.dat", "w") as f:
                f.write("#plumed.dat\n" +
                        "puck: PUCKERING ATOMS=" + 
                        c4 + ',' +
                        o4 + ',' +
                        c1 + ',' +
                        c2 + ',' +
                        c3 + "\n" + 
                        "PRINT ARG=puck.* FILE=" + colvar_name)

        except Exception as e:
            with open('/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/logs/log_de_erros.txt', 'a') as f:
                f.write(f"Erro: {str(e)}\n" + " file: " + ring_type + " - " + atoms_dict)

    return colvar_name

def run_puck(colvar_name):
    parts = colvar_name.split("/colvar/", 1)
    pdb_file_name = parts[1] if len(parts) > 1 else ""
    #Comando a ser executado
    command = "plumed driver --plumed puck.dat --mf_pdb ../pdb_sugars/" + pdb_file_name
    # Executar o comando e capturar a saída
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    print(result)

def separate_sugars(fileName):

    print("Separating sugars: " + fileName)
    
    #Cria um dicionário a partir do arquivo .cif
    mmcif_dict = MMCIF2Dict(fileName)

    try:
        #Separa as informações dos átomos
        atom_dict = {"group": mmcif_dict['_atom_site.group_PDB'], 
                        "id": mmcif_dict['_atom_site.id'], 
                        "label_atom_id":  mmcif_dict['_atom_site.label_atom_id'],
                        "label_comp_id":  mmcif_dict['_atom_site.label_comp_id'], 
                        "label_asym_id":  mmcif_dict['_atom_site.label_asym_id'],
                        "label_entity_id":  mmcif_dict['_atom_site.label_entity_id'],
                        "Cartn_x":  mmcif_dict['_atom_site.Cartn_x'],
                        "Cartn_y":  mmcif_dict['_atom_site.Cartn_y'],
                        "Cartn_z":  mmcif_dict['_atom_site.Cartn_z'],
                        "occupancy": mmcif_dict['_atom_site.occupancy'],
                        "B_iso_or_equiv": mmcif_dict['_atom_site.B_iso_or_equiv'],
                        "type_symbol": mmcif_dict['_atom_site.type_symbol']}
        #Cria o dataframe de átomos
        atom_df = pd.DataFrame(data = atom_dict)

        #Separa os hetero atomos
        hetatm_df = atom_df.loc[atom_df['group'] == 'HETATM']

        #Transforma as colunas numéricas em float
        hetatm_df["Cartn_x"] = hetatm_df["Cartn_x"].astype(float)
        hetatm_df["Cartn_y"] = hetatm_df["Cartn_y"].astype(float)
        hetatm_df["Cartn_z"] = hetatm_df["Cartn_z"].astype(float)
        hetatm_df["occupancy"] = hetatm_df["occupancy"].astype(float)
        hetatm_df["B_iso_or_equiv"] = hetatm_df["B_iso_or_equiv"].astype(float)

        #Separa so os carboidratos
        carbo_dict = pd.read_csv("/home/douglas/carboanalysis/carboanalysis/pdb/dicts/CCD_carbohydrate_list.tsv", sep = "\t", header = None, names = ['carbo_id', 'release_status'])
        carbo_list = carbo_dict["carbo_id"].values
        hetatm_df = hetatm_df.loc[hetatm_df['label_comp_id'].isin(carbo_list)]

        #reseta os index
        hetatm_df = hetatm_df.reset_index()

        #Itera os açúcares e escreve um arquivo no formato pdb pra cada um deles
        for index, row in hetatm_df.iterrows():

            if(index == 0):
                iter_sugar = row["label_comp_id"]
                iter_first_atom_id = row['id']
                sugar_atom_dict = {row["label_atom_id"]: row['id']}
                print(row["label_atom_id"] + " => " + row['id'])
                write_sugar_line_pdb(fileName, row, iter_sugar, iter_first_atom_id)

            else:
                #Verifica se é a prmeira linha
                if(index == len(hetatm_df.index) - 1):
                    write_sugar_line_pdb(fileName, row, iter_sugar, iter_first_atom_id)
                    #criar .dat e run bash
                    ring_type = is_piranose_or_furanose(iter_sugar, mmcif_dict)
                    colvar_name = alter_dat(sugar_atom_dict, ring_type, fileName, iter_sugar, iter_first_atom_id)
                    #roda plumed driver com o dat criado
                    run_puck(colvar_name)

                #Verifica se é o próximo açúcar
                if(row["label_comp_id"] != iter_sugar or row["label_atom_id"] == 'C1'):
                    #criar .dat e run bash
                    ring_type = is_piranose_or_furanose(iter_sugar, mmcif_dict)
                    colvar_name = alter_dat(sugar_atom_dict, ring_type, fileName, iter_sugar, iter_first_atom_id)
                    #roda plumed driver com o dat criado
                    run_puck(colvar_name)

                    sugar_atom_dict = {}
                    iter_sugar = row["label_comp_id"]
                    iter_first_atom_id = row["id"]
                
                write_sugar_line_pdb(fileName, row, iter_sugar, iter_first_atom_id)
                sugar_atom_dict[row["label_atom_id"]] = row['id']
    except:
        pass

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
        # print("Bfactoring...")
        # results = [i for i in pool.map(bfactor_values, fileNames) if i is not None]
        # print("Bfactor parsing Done!")
        print("Linkage...")
        results = [i for i in pool.map(find_linkages, fileNames) if i is not None]
        print("Linkage parsing Done!")

    print("thread time: ", time.time() - start)
