# HG_search
TATA.csv is an example file with PDB name and base pair which show features of mismodeled Watson-Crick.
1. pdb_rm_refine.sh is dowloading .pdb and .mtz listed in the .csv in folder PDB. Next, it chooses a purine in the base pair, delets it and copies in PDB_without_nt together with .mtz and refines it.
2. B_factor.py calculated mean B-factor for the nt of interest. Classify.sh moves structures to Bad_data if B-factor of the nt is >75.
(This will be changed to EDIA).

3. Next step is to place nt in WC and HG conformations, refine them and compare. 
