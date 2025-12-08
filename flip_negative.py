import sys
import glob
import os
import shutil
import pandas
from pymol import cmd


print("Flip A/G:", sys.argv[1],sys.argv[2],sys.argv[3])
pdb_name = str(sys.argv[1])
chain = str(sys.argv[2])
nt = str(sys.argv[3])
pdb_path = pdb_name
pdb_file = str(pdb_name + '.pdb')
cmd.load(pdb_file)
#print(str(pdb_name + "//" + chain + "/" + nt + "/O4'"),str(pdb_name + "//" + chain + "/" + "DC`" + nt + "/C1'"), str(pdb_name + "//" + chain + "/" + nt + "/N9"),str(pdb_name + "//" + chain + "/" + nt + "/C4"))
#x=cmd.get_dihedral("3pzp//Q/4/O4'","3pzp//Q/4/C1'","3pzp//Q/4/N9","3pzp//Q/4/C4")
a = f"{pdb_name}//{chain}/`\{nt}/O4'"

b = f"{pdb_name}//{chain}/`\{nt}/C1'"
c = f"{pdb_name}//{chain}/`\{nt}/N9"
d = f"{pdb_name}//{chain}/`\{nt}/C4"

x = cmd.get_dihedral(a, b, c, d)
print(x)
x += 180
cmd.set_dihedral(a, b, c, d, x)

cmd.save(pdb_name + '_flipped.pdb', pdb_name)
