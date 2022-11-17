import os
import pandas as pd
import warnings
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

    for fileName in fileNames:
        
        mmcif_dict = MMCIF2Dict(fileName)
        chem_components = mmcif_dict["_chem_comp.id"]

        if(set(carbo_dict["carbo_id"].values) & set(chem_components)):
            
            filteredFileNames.append(fileName)
            filteredEntries.append(mmcif_dict["_entry.id"][0])

    # escreve um arquivo .txt com o nome das entradas cujas estrutuaras possuam carbohidratos
    with open("/home/douglas_lima/pdb/dataframes/carbo_entrys.txt", "w") as file:
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
        entity_number_of_molecules = mmcif_dict["_entity.pdbx_number_of_molecules"]
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

            carbo_dict = pd.read_csv("/home/douglas_lima/pdb/dicts/CCD_carbohydrate_list.tsv", sep = "\t", header = None, names = ['carbo_id', 'REF'])

            monosaccharide_df = monosaccharide_df[monosaccharide_df.comp_id.isin(carbo_dict["carbo_id"].values)]
           
        # junta os dataframes dos monossacarídeos
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
    
    monosaccharides.to_csv(path_or_buf="/home/douglas_lima/pdb/dataframes/monosaccharides.csv")
    oligosaccharide_df.to_csv(path_or_buf="/home/douglas_lima/pdb/dataframes/oligosaccharides.csv")

os.chdir("/home/douglas_lima/pdb/testesCif")
fileNames = os.listdir("/home/douglas_lima/pdb/testesCif")

filtrados = filter_maxResolution(fileNames, 1.6)
filtrados = filter_maxOWAB(fileNames, 60)
filtrados = filter_structureMethod(fileNames, 'X-RAY CRYSTALLOGRAPHY')
separate(fileNames)
