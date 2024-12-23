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
from tqdm import tqdm

def is_piranose_or_furanose(sugar_id, mmcif_dict):
    try:
        chem_comp_dict = {"comp_id": mmcif_dict['_chem_comp.id'], "name": mmcif_dict['_chem_comp.name']}
        chem_comp_df = pd.DataFrame(data = chem_comp_dict)

        if("pyranose" in chem_comp_df[chem_comp_df.comp_id == sugar_id]["name"].values[0]):
            #print("pyr")
            return "pyranose"

        elif("furanose" in chem_comp_df[chem_comp_df.comp_id == sugar_id]["name"].values[0]):
            #print("fur")
            return "furanose"
    except Exception as e:
        with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/logs/puckering_log.txt", "a") as f:
            f.write("Exception pir fur in : " + mmcif_dict['_entry.id'] + " =>" + str(e) + "\n")
        return None

def write_sugar_line_pdb(fileName, row, sugar, sugar_first_id, new_atom_id):
    try:
        with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/missing_sugars/" + fileName[:4] + "_" + sugar + "_" + sugar_first_id + ".pdb", "a") as f:
            # Formatar a linha de acordo com o padrão PDB
            pdb_line = "{:<6}{:>5} {:<4}{:<3} {:<2}{:>4}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:>2}".format(
                row['group'],            # 'HETATM'
                new_atom_id,             # Número serial do átomo
                row['label_atom_id'],    # Nome do átomo
                row['label_comp_id'],    # Nome do resíduo
                row['label_asym_id'],    # ID da cadeia
                row['label_entity_id'],  # Número da sequência do resíduo
                row['Cartn_x'],          # Coordenada x
                row['Cartn_y'],          # Coordenada y
                row['Cartn_z'],          # Coordenada z
                row['occupancy'],        # Ocupação
                row['B_iso_or_equiv'],   # Fator de temperatura
                row['type_symbol']       # Elemento químico
            )

            f.write(pdb_line + '\n')
    except Exception as e:
        with open('/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/logs/log_de_erros.txt', 'a') as f:
            f.write("Exception in " + fileName + ": " + str(e) + "\n")
    
        return None

#Adiciona comentário no PDB
def add_remarks(ring_atoms, fileName, sugar, sugar_first_id):
    with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/missing_sugars/" + fileName[:4] + "_" + sugar + "_" + sugar_first_id + ".pdb", "a") as file:
        file.write(f"REMARK 1 {ring_atoms}\n")
#Retorna a lista de átomos do anel
def ring_atoms_identifier(atoms_dict, ring_type, fileName):
    if(ring_type == "pyranose"):
        try:
            o5 = atoms_dict['O5']
            c1 = atoms_dict['C1']
            c2 = atoms_dict['C2']
            c3 = atoms_dict['C3']
            c4 = atoms_dict['C4']
            c5 = atoms_dict['C5']
            #print(atoms_dict)
            return o5 + "," + c1 + "," + c2 + "," + c3 + "," + c4 + "," + c5
        
        except Exception as e:
            with open('/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/logs/log_de_erros_ht2.txt', 'a') as f:
                f.write("Exception: pyranose => " + fileName +  "\n")
            return None

    elif(ring_type == "furanose"):
        try:
            c4 = atoms_dict['C4']
            o4 = atoms_dict['O4']
            c1 = atoms_dict['C1']
            c2 = atoms_dict['C2']
            c3 = atoms_dict['C3']
            #print(atoms_dict)
            return c4 + "," + o4 + "," + c1 + "," + c2 + "," + c3
        
        except Exception as e:
            with open('/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/logs/log_de_erros_ht2.txt', 'a') as f:
                f.write("Exception: pyranose => " + fileName +  "\n")
            return None
    else:
        #print('nada')
        pass

def read_remmarks(fileName, remark_number):
    comments = []
    remark_prefix = f"REMARK {remark_number}"
    with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/missing_sugars/" + fileName, 'r') as file:
        for line in file:
            if line.startswith(remark_prefix):
                comments.append(line.strip())
    return comments

def parse_pdb_line(line):
    return {
        'group': line[0:6].strip(),
        'atom_id': int(line[6:11].strip()),
        'label_atom_id': line[12:16].strip(),
        'label_comp_id': line[17:20].strip(),
        'label_asym_id': line[21].strip(),
        'label_entity_id': line[22:26].strip(),
        'Cartn_x': float(line[30:38].strip()),
        'Cartn_y': float(line[38:46].strip()),
        'Cartn_z': float(line[46:54].strip()),
        'occupancy': float(line[54:60].strip()),
        'B_iso_or_equiv': float(line[60:66].strip()),
        'type_symbol': line[76:78].strip()
    }

def read_pdb_to_dataframe(file_path):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith(('ATOM', 'HETATM')):
                parsed_line = parse_pdb_line(line)
                data.append(parsed_line)
    return pd.DataFrame(data)

def puck_subprocess(fileName):
    try:
        #Comando a ser executado
        command = "plumed driver --plumed puck.dat --mf_pdb ../missing_sugars/" + fileName
        # Executar o comando e capturar a saída
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
    except Exception as e:
        with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/logs/puckering_missing_log.txt", "a") as f:
            f.write("PUCK SUBPROCESS Exception in " + fileName + ": " + str(e) + "\n")
        return None

def puck_calcs(fileName):
    try:
        remmarks = read_remmarks(fileName, 1)
        remmark = remmarks[0]
        ring_atoms = remmark[9:]
        with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/puckering/puck.dat", "w") as f:
            f.write("#plumed.dat \n" +
                    "puck: PUCKERING ATOMS="+ ring_atoms + "\n" + 
                    "PRINT ARG=puck.* FILE=/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/puckering/colvar_missing_sugars/" + fileName[:-4])
        
        puck_subprocess(fileName)
    except Exception as e:
        with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/logs/only_calc_puck_missing_errors.txt", "a") as f:
                    f.write("PUCK_CALCS Error in "+ fileName +": " + e)


def only_separate_sugars(fileName):
    try:
        os.chdir("/home/douglas/carboanalysis/data/unzipped")
        # Verifica se o arquivo existe
        if not os.path.isfile(fileName):
            raise FileNotFoundError(f"File not found: {fileName}")

        #print("Separating sugars: " + fileName)
        
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
                            "type_symbol": mmcif_dict['_atom_site.type_symbol'],
                            "auth_seq_id": mmcif_dict['_atom_site.auth_seq_id'],
                            "label_asym_id": mmcif_dict['_atom_site.label_asym_id'],
                            "auth_asym_id": mmcif_dict['_atom_site.auth_asym_id']}
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
                    iter_auth_seq_id = row['auth_seq_id']
                    iter_label_asym_id = row['label_asym_id']
                    iter_auth_asym_id = row['auth_asym_id']
                    atom_id_cont = 1
                    sugar_atom_dict = {row["label_atom_id"]: str(atom_id_cont)}

                    write_sugar_line_pdb(fileName, row, iter_sugar, iter_first_atom_id, str(atom_id_cont))

                else:
                    #Verifica se é a última linha
                    if(index == len(hetatm_df.index) - 1):
                        atom_id_cont += 1
                        write_sugar_line_pdb(fileName, row, iter_sugar, iter_first_atom_id, str(atom_id_cont))
                        sugar_atom_dict[row["label_atom_id"]] = str(atom_id_cont)
                        #Separa e adiciona comentário
                        ring_type = is_piranose_or_furanose(iter_sugar, mmcif_dict)
                        
                        ring_atoms = ring_atoms_identifier(sugar_atom_dict, ring_type, fileName)

                        add_remarks(ring_atoms, fileName, iter_sugar, iter_first_atom_id)
                        
                        return 0

                    #Verifica se é o próximo açúcar
                    #if(row["label_comp_id"] != iter_sugar or row["label_atom_id"] == 'C1'):
                    if((row["label_comp_id"] != iter_sugar) or (row["auth_seq_id"] != iter_auth_seq_id) or (row['label_asym_id'] != iter_label_asym_id) or (row['auth_asym_id'] != iter_auth_asym_id)):
                        #Separa e adiciona comentário
                        ring_type = is_piranose_or_furanose(iter_sugar, mmcif_dict)
                        
                        ring_atoms = ring_atoms_identifier(sugar_atom_dict, ring_type, fileName)

                        add_remarks(ring_atoms, fileName, iter_sugar, iter_first_atom_id)

                        sugar_atom_dict = {}
                        iter_sugar = row["label_comp_id"]
                        iter_first_atom_id = row["id"]
                        iter_auth_seq_id = row['auth_seq_id']
                        iter_label_asym_id = row['label_asym_id']
                        iter_auth_asym_id = row['auth_asym_id']
                        atom_id_cont = 0
                    
                    atom_id_cont += 1
                    write_sugar_line_pdb(fileName, row, iter_sugar, iter_first_atom_id, str(atom_id_cont))
                    sugar_atom_dict[row["label_atom_id"]] = str(atom_id_cont)
        except Exception as e:
            # Registra o erro em um log específico para esta função
            with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/logs/separate_sugars_errors_only_missing.txt", "a") as f:
                f.write("Exception in " + fileName + ": " + str(e) + "\n")
            return None
    except Exception as e:
        # Registra o erro em um log específico para esta função
        with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/logs/separate_sugars_errors_only_missing.txt", "a") as f:
            f.write("Exception in " + fileName + ": " + str(e) + "\n")
        return None

def check_remmarks(fileName):
    remmark = read_remmarks(fileName, 1)
    if remmark:
        if 'None' in remmark[0]:
            print(f'{fileName} -> {remmark}')
            return fileName
        
def update_remmarks(fileName):
    pass

if __name__ == "__main__":

    # df = pd.read_csv("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/missing_entries.csv", sep=';')
    # df["entry_id"] = df['entry_id'].str.lower()
    # df["entry_id"] = df['entry_id'] + '.cif'
    # fileNames = df['entry_id'].values
    
    # os.chdir("/home/douglas/carboanalysis/data/unzipped")

    os.chdir("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/missing_sugars/")
    fileNames = os.listdir('/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/missing_sugars/')
    with Pool() as pool:
        # print("Separating sugars...")
        # results = []
        # for result in tqdm(pool.imap(only_separate_sugars, fileNames), total=len(fileNames)):
        #     if result is not None:
        #         results.append(result)
        # print("Separating Done!")
        print("Reading remmarks...")
        results = []
        for result in tqdm(pool.imap(check_remmarks, fileNames), total=len(fileNames)):
            if result is not None:
                results.append(result)
        print("Done!")
        

    # file_names = os.listdir("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/pdb_res_higher_than_2_sugars/")
    # os.chdir("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/puckering")
    # for file_name in tqdm(file_names, total=len(file_names)):
    #     puck_calcs(file_name)