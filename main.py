import Bio
import numpy as np
from Bio.PDB import PDBList
import matplotlib.pyplot as plt
import sys
import warnings

pdb_list = PDBList()

arg1 = sys.argv[1] #nazwa bialka
arg2 = sys.argv[2] #ktory lancuch
arg3 = sys.argv[3] #odleglosc

#tutaj sb pobieram wpisane w konsoli białko
pdb_id = str(arg1)

warnings.simplefilter("ignore")
pdb_filename = pdb_list.retrieve_pdb_file(pdb_id, pdir="data/PDB_files", file_format="pdb")
structure = Bio.PDB.PDBParser().get_structure(pdb_id, pdb_filename)
warnings.resetwarnings()

atom_coords = []

chain = structure[0][arg2]

chain.atom_to_internal_coordinates(verbose=True)
counter = 0
#do atom_array, trafiaja tylko i wylacznie koordynaty atomow o identyfikatorze CA
for res in chain:
    for atom in res:
        if atom.get_id() == "CA":
            counter+=1
            atom_coords.append(atom.get_coord())

#zmieniam na numpy array
atom_coords = np.array(atom_coords)

#tu robie distance matrix wypelniony zerami, zeby go potem zapelnic
matrix_bool = np.zeros((len(atom_coords), len(atom_coords)), dtype=bool)

for i in range(len(atom_coords)):
    for j in range(len(atom_coords)):
        distance = np.linalg.norm(atom_coords[i] - atom_coords[j])
        #8 oznacza dystans ktory jest nasza granica - potem to trzeba zmienic na wpis z konsoli
        if distance < int(arg3):
            matrix_bool[i, j] = True
        else:
            matrix_bool[i, j] = False

#teraz trzeba zrobic distances i w zalenosci od okreslonej odleglosci zrobic macierz boolowska i wyprintowac plot

plt.imshow(matrix_bool, cmap="gray")
plt.show()














