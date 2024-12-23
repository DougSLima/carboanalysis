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
            entity_dict = {"entity_id": mmcif_dict['_entity.id'], 
                    "type": mmcif_dict['_entity.type'], 
                   "description": mmcif_dict['_entity.pdbx_description'], 
                    "n_of_molecules": mmcif_dict['_entity.pdbx_number_of_molecules']}

            entity_df = pd.DataFrame(data=entity_dict)
            entity_df['entry_id'] = entry_id
            entity_df = entity_df.loc[entity_df['type'] == 'branched']

            if entity_df.empty:
                pass
            else: 
                entity_df.to_csv("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/final_analysis/branched.csv", header=None, sep=';', mode='a')

            return entity_df
        except Exception as e:
            # Registra o erro em um log específico para esta função
            with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/logs/get_branched_log.txt", "a") as f:
                f.write("Exception in " + entry_id + " - res: " + str(e) + "\n")
            return None

def get_oligos(entry_id):

    mmcif_dict = MMCIF2Dict(entry_id + ".cif")

    if '_pdbx_entity_branch_list.entity_id' in mmcif_dict.keys():
        try:
            oligo_dict = {"entity_id": mmcif_dict["_pdbx_entity_branch_list.entity_id"], 
                    "comp_id": mmcif_dict["_pdbx_entity_branch_list.comp_id"], 
                   "num": mmcif_dict["_pdbx_entity_branch_list.num"]}

            oligo_df = pd.DataFrame(data=oligo_dict)

            oligo_df['entry_id'] = entry_id

            oligo_df['num'] = pd.to_numeric(oligo_df['num'])

            result_df = oligo_df.groupby('entity_id', as_index=False)['num'].max()

            branched_df = get_branched(entry_id)
            
            if branched_df.empty:
                return None

            branched_df = branched_df[['entity_id', 'n_of_molecules']]

            result_df = result_df.merge(branched_df, how='left', left_on='entity_id', right_on='entity_id')
            result_df['entry_id'] = entry_id 
            result_df = result_df[['entry_id', 'entity_id', 'num', 'n_of_molecules']]
            
            return result_df
        except Exception as e:
            # Registra o erro em um log específico para esta função
            with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/logs/get_oligo_log.txt", "a") as f:
                f.write("Exception in " + entry_id + " - res: " + str(e) + "\n")
            return None

def write_oligos(entry_id):
    oligo_df = get_oligos(entry_id)
    if oligo_df is not None:
        oligo_df.to_csv("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/final_analysis/oligos.csv", header=None, sep=';', mode='a')

def get_ramifications():

    names = ['entry_id', 'entity_id', 'num', 'n_of_molecules']
    oligo_df = pd.read_csv("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/final_analysis/oligos.csv", header=None, sep=';', names = names)

    for index, row in oligo_df.iterrows():
        try: 
            mmcif_dict = MMCIF2Dict(row['entry_id'] + ".cif")

            if '_pdbx_entity_branch_link.link_id' in mmcif_dict.keys():
                branch_link_dict = {"entity_id": mmcif_dict["_pdbx_entity_branch_link.entity_id"], 
                    "num_1": mmcif_dict["_pdbx_entity_branch_link.entity_branch_list_num_1"], 
                   "comp_id_1": mmcif_dict["_pdbx_entity_branch_link.comp_id_1"],
                   "num_2": mmcif_dict["_pdbx_entity_branch_link.entity_branch_list_num_2"], 
                   "comp_id_2": mmcif_dict["_pdbx_entity_branch_link.comp_id_2"]}

                branch_link_df = pd.DataFrame(data=branch_link_dict)
                #Seleciona somente a entity_id de entrada
                branch_link_df['entity_id'] = branch_link_df['entity_id'].astype(int)

                branch_link_df = branch_link_df.loc[branch_link_df['entity_id'] == row['entity_id']]
                
                branch_link_df['col1'] = branch_link_df['num_1'] + branch_link_df['comp_id_1']  
                branch_link_df['col2'] = branch_link_df['num_2'] + branch_link_df['comp_id_2']
                
                comps = list(branch_link_df['col1'].values) + list(branch_link_df['col2'].values)
                comps = pd.Series(comps)
                comps = comps.value_counts()
                count_gte_3 = (comps >= 3).sum()
                if(count_gte_3 > 0):
                    ram_df = pd.DataFrame([[row['entry_id'], row['entity_id'], row['num'], row['n_of_molecules'], True, count_gte_3]],
                    columns=['entry_id', 'entity_id', 'num', 'n_of_molecules', 'ramification', 'n_of_ramifications'])
                    ram_df.to_csv("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/final_analysis/ramifications.csv", header=None, sep=';', mode='a')
                else:
                    ram_df = pd.DataFrame([[row['entry_id'], row['entity_id'], row['num'], row['n_of_molecules'], False, 0]],
                    columns=['entry_id', 'entity_id', 'num', 'n_of_molecules', 'ramification', 'n_of_ramifications'])
                    ram_df.to_csv("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/final_analysis/ramifications.csv", header=None, sep=';', mode='a')

        except Exception as e:
            # Registra o erro em um log específico para esta função
            with open("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/logs/get_ramifications_log.txt", "a") as f:
                f.write("Exception in " + row['entry_id'] + " - res: " + str(e) + "\n")
            return None


if __name__ == "__main__":

    df = pd.read_csv('/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/carbo_entrys_nd.txt', names = ['entry_filename'])
    fileNames = df['entry_filename'].str.replace('.cif', '').values

    os.chdir("/home/douglas/carboanalysis/data/unzipped")

    # with Pool() as pool:
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
        # print("Organism and PolyChainSize info...")
        # results = []
        # for result in tqdm(pool.imap(write_org_polysize, fileNames), total=len(fileNames)):
        #     if result is not None:
        #         results.append(result)
        # print("Done!")
        # print("Get oligo...")
        # results = []
        # for result in tqdm(pool.imap(write_oligos, fileNames), total=len(fileNames)):
        #     if result is not None:
        #         results.append(result)
        # print("Done!")

    get_ramifications()