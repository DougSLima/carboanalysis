import os
from Bio.PDB import *

pdbParser = PDBParser()

os.chdir("/home/douglas_lima/pdb/testes")
fileNames = os.listdir("/home/douglas_lima/pdb/testes")

for fileName in fileNames:
    file = open(fileName, 'r')
    structure = pdbParser.get_structure(file.name, file)
    print(structure.header.get("resolution"))

# for fileName in fileNames:
#     file = open(fileName, 'r')
#     structure = pdbParser.get_structure(file.name, file)
#     atomList = Selection.unfold_entities(structure, 'A')
#     occupancy_x_bfactor = 0

#     for atom in atomList:
#         occupancy_x_bfactor += (atom.get_occupancy() * atom.get_bfactor())
    
#     owab = occupancy_x_bfactor/len(atomList)
#     print("Nome: " + str(fileName) + " - OWAB: " + str(owab))






