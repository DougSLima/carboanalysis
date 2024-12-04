import pandas as pd

from SIC2023_resolution import *

def ion_find(fileName):
    
    print("Ion finding...: " + fileName)

    mmcif_dict = MMCIF2Dict(fileName)

    entity_dict = {"id": mmcif_dict['_entity.id'], 
                   "type": mmcif_dict['_entity.type'], 
                   "description": mmcif_dict['_entity.pdbx_description'], 
                    "number_of_molecules": mmcif_dict['_entity.pdbx_number_of_molecules']}
    
    entity_df = pd.DataFrame(data = entity_dict)

    entity_df = entity_df[entity_df['type'] == 'non-polymer']

    entity_df = entity_df[entity_df['description'].str.contains('ION', na=False)]

    if not entity_df.empty:
        ions = dict(zip(entity_df['description'].tolist(), entity_df['number_of_molecules'].tolist()))
        
        ion_dict = {
            'entry': fileName,
            'ion': ions
        }

        df = pd.DataFrame([ion_dict])  # Crie um DataFrame com uma lista de um único dicionário
        df.to_csv('/home/douglas/ions.csv', mode='a', header=False, index=False, sep=';')

if __name__ == '__main__':

    df = pd.read_csv('/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/puckering/piranoses_clean_family_and_protein_name_entry_name.csv')
    df = df.iloc[:, 4:]
    # Remover todos os números da coluna
    df['Family'] = df['Family'].str.replace('\d+', '', regex=True)

    cbm_df = df[df['Family'] ==  'CBM']

    cbm_entries = (df['sugar_entry'] + '.cif').drop_duplicates()

    file_names = cbm_entries.values

    os.chdir("/home/douglas/carboanalysis/data/unzipped")

    with Pool() as pool:
        print("Ion finding...")
        results = [i for i in pool.map(ion_find, file_names) if i is not None]
        print("Done!")

