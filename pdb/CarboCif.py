import os
from Bio.PDB import *
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning #ignorar warning (PDBConstructionWarning: WARNING: Chain B is discontinuous at line numeroDaLinha.)
from collections import defaultdict

warnings.simplefilter('ignore', PDBConstructionWarning)

import os
from Bio.PDB import *
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning #ignorar warning (PDBConstructionWarning: WARNING: Chain B is discontinuous at line numeroDaLinha.)
from collections import defaultdict

warnings.simplefilter('ignore', PDBConstructionWarning)

# Remove da lista recebida todos os nomes de arquivos cujas estruturas apresentem resoluções maiores que a resolução máxima selecionada
def filter_maxResolution(fileNames, maxResolution):

    filteredFileNames = [] # lista de filtrados
    pdbParser = MMCIFParser()

    for fileName in fileNames:
        file = open(fileName, 'r')
        structure = pdbParser.get_structure(file.name, file)
        
        if structure.header.get("resolution") <= maxResolution: # Testa se a resolução da estrutura é menor ou igual a resolução desejada (atributo da função)
            filteredFileNames.append(fileName)
    
    return filteredFileNames

# Remove da lista recebida todos os nomes de arquivos cujas estruturas apresentem OWAB maiores que o OWAB máximo selecionado
def filter_maxOWAB(fileNames, maxOWAB):

    filteredFileNames = []
    pdbParser = MMCIFParser()

    for fileName in fileNames:
        file = open(fileName, 'r')
        structure = pdbParser.get_structure(file.name, file)
        atomList = Selection.unfold_entities(structure, 'A')
        occupancy_x_bfactor = 0

        for atom in atomList:
            occupancy_x_bfactor += (atom.get_occupancy() * atom.get_bfactor())
        
        owab = occupancy_x_bfactor/len(atomList) # OWAB = média de ocupância * b_factor dos átomos da estrutura
        print(fileName)
        print(owab)
        if owab <= maxOWAB:
            filteredFileNames.append(fileName)
    
    return filteredFileNames




os.chdir("/home/douglas_lima/pdb/testesCif")
fileNames = os.listdir("/home/douglas_lima/pdb/testesCif")

print(fileNames)
print(filter_maxResolution(fileNames, 3))