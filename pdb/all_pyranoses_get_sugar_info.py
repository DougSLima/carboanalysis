#Get info
from SIC2023_resolution import *
import numpy as np
from tqdm import tqdm
import pandas as pd
from multiprocessing import Pool, Lock
from tqdm import tqdm

##Piranoses get info
# xray_filtered_piranose_df = pd.read_csv('/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/carbovault/all_monossacharide_df_without_missing.csv', header=None, sep=';',
#                                         names = ['sugar','iupac_name', 'puck.phi_graus', 'puck.theta_graus', 'entry_id','entry_resolution', 'is_in_cazy', 'type'])
# em_piranose_df = pd.read_csv("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/puckering/EM_piranoses_degrees.csv")
# ht2_piranose_df = pd.read_csv("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/puckering/ht2_piranoses_degrees.csv", header=None,
#                         names = ['sugar', 'time', 'puck.qx', 'puck.qy', 'puck.qz', 'puck.phi', 'puck.theta', 'puck.amplitude', 'puck.phi_graus', 'puck.theta_graus'])

# #Remove unnamed columns
# em_piranose_df = em_piranose_df.iloc[:, 1:]

# #pega os sugar ids
# xray_filtered_piranoses = xray_filtered_piranose_df['sugar'].values
# em_piranoses = em_piranose_df['sugar'].values
# ht2_piranoses = ht2_piranose_df['sugar'].values

# #concatena
# all_piranoses = np.concatenate([xray_filtered_piranoses, em_piranoses, ht2_piranoses])

##Furanoses Get info

xray_filtered_furanose_df = pd.read_csv('/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/puckering/furanoses_antesdebfactor_mean.csv')
xray_filtered_furanose_df['sugar'] = xray_filtered_furanose_df['sugar'].str.split('filtered/').str[-1]
em_furanose_df = pd.read_csv("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/puckering/EM_furanoses.csv")
ht2_furanose_df = pd.read_csv("/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/puckering/ht2_furanoses.csv")

all_furanoses = list(xray_filtered_furanose_df['sugar'].values) + list(em_furanose_df['sugar'].values) + list(ht2_furanose_df['sugar'].values)

os.chdir("/home/douglas/carboanalysis/data/unzipped")

# Função wrapper para incluir o output_path
def get_sugar_info_with_output_path(sugar_id):
    #return get_sugar_info(sugar_id, '/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/puckering/all_pyranoses/all_pyranoses.csv')
    return get_sugar_info(sugar_id, '/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/puckering/all_pyranoses/all_furanoses.csv')

with Pool() as pool:
    print("Get sugar info...")
    results = []
    for result in tqdm(pool.imap(get_sugar_info_with_output_path, all_furanoses), total=len(all_furanoses)):
        if result is not None:
            results.append(result)
    print("Done!")
