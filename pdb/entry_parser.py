from multiprocessing import Process, Pool
from tqdm import tqdm
import pandas as pd
import logging
from SIC2023_resolution import * 

# Configuração básica do log
logging.basicConfig(
    filename='/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/logs/entry_info.log',
    level=logging.ERROR,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def parse_entry(file_path):
    
    try:
        mmcif_dict = MMCIF2Dict(file_path)

        # Pega as informações das entradas
        entry_dict = {"entry_id": mmcif_dict['_entry.id'], 
                        "initial_deposition_date": mmcif_dict['_pdbx_database_status.recvd_initial_deposition_date'],
                        "status_code": mmcif_dict['_pdbx_database_status.status_code'], 
                        "method":  mmcif_dict['_exptl.method'], 
                        "title": mmcif_dict["_struct.title"], 
                        "keywords": mmcif_dict['_struct_keywords.pdbx_keywords'], 
                        "other_keywords": mmcif_dict['_struct_keywords.text']}
        entry_df = pd.DataFrame(data = entry_dict)

        entry_df.to_csv('/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/entry_info.csv', mode='a', header=None, sep=';')
    except Exception as e:
        logging.error(f"Erro no PDBID {file_path}: {e}")

def write_resolution(file_name: str):
    res = get_resolution_all(file_name)
    if res:
        name = file_name.replace('.cif', '')
        name = file_name.upper()
        
        res_df = pd.DataFrame([[name, res]], columns=["entry_id", "resolution"])
        res_df.to_csv(path_or_buf="/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/resolution_all.csv", mode='a', index=False, header=False, sep=';')

def get_resolution(entry_id):
    
    mmcif_dict = MMCIF2Dict(entry_id + ".cif")

    if '_refine.ls_d_res_high' in mmcif_dict.keys():
        try:
            return mmcif_dict["_refine.ls_d_res_high"][0]
        except Exception as e:
            # Registra o erro em um log específico para esta função
            with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/logs/get_resolution_log.txt", "a") as f:
                f.write("Exception in " + entry_id + " - res: " + mmcif_dict["_refine.ls_d_res_high"][0] + ": " + str(e) + "\n")
            return None


def get_organism(entry_id):

    mmcif_dict = MMCIF2Dict(entry_id + ".cif")

    if '_entity_src_gen.pdbx_gene_src_scientific_name' in mmcif_dict.keys():
        try:
            return mmcif_dict["_entity_src_gen.pdbx_gene_src_scientific_name"][0]
        except Exception as e:
            # Registra o erro em um log específico para esta função
            with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/logs/get_organism_log.txt", "a") as f:
                f.write("Exception in " + entry_id + " - res: " + mmcif_dict["_entity_src_gen.pdbx_gene_src_scientific_name"][0] + ": " + str(e) + "\n")
            return None
        
def get_poly_length(entry_id):

    mmcif_dict = MMCIF2Dict(entry_id + ".cif")

    if '_entity_poly_seq.entity_id' in mmcif_dict.keys():
        try:
            return len(mmcif_dict["_entity_poly_seq.entity_id"])
        except Exception as e:
            # Registra o erro em um log específico para esta função
            with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/logs/get_poly_length_log.txt", "a") as f:
                f.write("Exception in " + entry_id + " - res: " + mmcif_dict["_entity_poly_seq.entity_id"] + ": " + str(e) + "\n")
            return None

def write_org_polysize(entry_id):
    organism = get_organism(entry_id)
    poly_length = get_poly_length(entry_id)
    if organism and poly_length:
        name = entry_id.upper()
        
        res_df = pd.DataFrame([[name, organism, poly_length]], columns=["entry_id", "organism", "poly_length"])
        res_df.to_csv(path_or_buf="/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/organism_and_polylength.csv", mode='a', index=False, header=False, sep=';')

def get_branched(entry_id):

    mmcif_dict = MMCIF2Dict(entry_id + ".cif")

    if '_entity.type' in mmcif_dict.keys():
        try:
            

            return len(mmcif_dict["_entity_poly_seq.entity_id"])
        except Exception as e:
            # Registra o erro em um log específico para esta função
            with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/logs/get_poly_length_log.txt", "a") as f:
                f.write("Exception in " + entry_id + " - res: " + mmcif_dict["_entity_poly_seq.entity_id"] + ": " + str(e) + "\n")
            return None



if __name__ == "__main__":

    df = pd.read_csv('/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/carbo_entrys_nd.txt', names = ['entry_filename'])
    fileNames = df['entry_filename'].str.replace('.cif', '').values

    os.chdir("/home/douglas/carboanalysis/data/unzipped")

    with Pool() as pool:
        # print("Entry info...")
        # results = []
        # for result in tqdm(pool.imap(parse_entry, fileNames), total=len(fileNames)):
        #     if result is not None:
        #         results.append(result)
        # print("Entry info Done!")
        # print("Res info...")
        # results = []
        # for result in tqdm(pool.imap(write_resolution, fileNames), total=len(fileNames)):
        #     if result is not None:
        #         results.append(result)
        # print("Done!")
        print("Organism and PolyChainSize info...")
        results = []
        for result in tqdm(pool.imap(write_org_polysize, fileNames), total=len(fileNames)):
            if result is not None:
                results.append(result)
        print("Done!")

