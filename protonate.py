import sys
import glob
import os
import shutil
import pandas
from pymol import cmd


print("Protonate C N3:", sys.argv[1],sys.argv[2],sys.argv[3])
pdb_name = str(sys.argv[1])
chain_C = str(sys.argv[2])
nt = str(sys.argv[3])
pdb_path = pdb_name
pdb_file = str(pdb_name + '.pdb')

cmd.load(pdb_file)


CN3 = str(pdb_name + "//" + chain_C + "/" + nt + "/N3")

cmd.h_add(selection=CN3 )

cmd.save(pdb_name + '_protonated.pdb', pdb_name)
	## add command to move flipped pdb into the flipped directory
