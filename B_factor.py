from Bio.PDB import *
import numpy as np
import pandas as pd

data = pd.read_csv("TATA.csv", sep = ';', header=None)
pdb = list(data[0])
b_factors = []
head = []
for i in range(len(pdb)):
    parser = PDBParser();
    structure = parser.get_structure(pdb[i], "PDB/"+pdb[i]+".pdb")
    atoms = structure.get_atoms()
    if data[14][i] == "DC" or data[14][i] == "DT":
        c = data[13][i];
        r = data[11][i][6:]
        h = data[11][i]
        print(c,r)
    else:
        c =  data[12][i]
        r = data[10][i][6:]
        h = data[10][i]
        print(c,r)
    residue = structure[0][str(c)][int(r)]
    b =[]
    for atom in residue:
        b.append(atom.bfactor)
    b_factors.append(int(np.mean(b)))
    head.append(str(pdb[i])+'_'+str(h))
print(head)
    
f = open('b_factor.txt', 'w')
for i in range(len(pdb)):
    print( pdb[i], b_factors[i])
    f.writelines(str(head[i]))
    f.writelines('\t')
    f.writelines(str(b_factors[i]))
    f.writelines('\n')
f.close()