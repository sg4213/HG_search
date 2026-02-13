#!/bin/bash
#rm -r PDB_without_nt   source /Users/sima/phenix-1.20.1-4487/phenix_env.csh
mkdir PDB_without_nt

if [ $# -ne 9 ]; then
    echo "Usage: $0 xxxx xxxx chain_1 nt_type_1 nt_number_1 chain_2 nt_type_2 nt_number_2 xxxx"
fi

pdb_id="$1"
pdb_file="$1_final.pdb"
mtz_file="$2_final.mtz"
chain_1="$3"
nt_type_1="$4"
nt_number_1="$5"
chain_2="$6"
nt_type_2="$7"
nt_number_2="$8"

echo "PDB File: $pdb_id"
echo "PDB File: $pdb_file"
echo "MTZ File: $mtz_file"
echo "Chain 1: $chain_1"
echo "Nucleotide Type 1: $nt_type_1"
echo "Nucleotide Number 1: $nt_number_1"
echo "Chain 2: $chain_2"
echo "Nucleotide Type 2: $nt_type_2"
echo "Nucleotide Number 2: $nt_number_2"

select_nt() {
  if [ "$5" == "G" ]
  then
       echo "$4 $6" >> delete.txt
       chain_delete=$chain_2
       nt_delete=$nt_number_2
  elif [ "$2" == "G" ]
  then
       echo "$1 $3" >> delete.txt
       chain_delete=$chain_1
       nt_delete=$nt_number_1
  elif [ "$5" == "A" ]
  then
       echo "$4 $6" >> delete.txt
       chain_delete=$chain_2
       nt_delete=$nt_number_2
  elif [ "$2" == "A" ]
  then
       echo "$1 $3" >> delete.txt
       chain_delete=$chain_1
       nt_delete=$nt_number_1
  else    echo "mismatch or modification"
  fi
}
########################################
#Select nt and copy files for refinment
########################################

select_nt $chain_1 $nt_type_1 $nt_number_1 $chain_2 $nt_type_2 $nt_number_2
mkdir PDB_without_nt/$pdb_id\_$chain_delete\_$nt_delete/
mkdir PDB_without_nt/$pdb_id\_$chain_delete\_$nt_delete/HG
mkdir PDB_without_nt/$pdb_id\_$chain_delete\_$nt_delete/WC
cp $mtz_file PDB_without_nt/$pdb_id\_$chain_delete\_$nt_delete
cp $pdb_file PDB_without_nt/$pdb_id\_$chain_delete\_$nt_delete/HG
cp $pdb_file PDB_without_nt/$pdb_id\_$chain_delete\_$nt_delete/WC

if [ $nt_delete -lt 0 ]
then
    cp flip_negative.py PDB_without_nt/$pdb_id\_$chain_delete\_$nt_delete/HG
elif [ $nt_delete -ge 0 ]
then
    cp flip.py PDB_without_nt/$pdb_id\_$chain_delete\_$nt_delete/HG
else
    cp flip.py PDB_without_nt/$pdb_id\_$chain_delete\_$nt_delete/HG
fi


########################################
# Check nt occupancy
########################################
cp check_occupancy.py PDB_without_nt/$pdb_id\_$chain_delete\_$nt_delete/WC
cd PDB_without_nt/$pdb_id\_$chain_delete\_$nt_delete/WC
python3 check_occupancy.py ${pdb_id}_final.pdb $chain_delete $nt_delete >> occupancy

occ="$(tr -d '[:space:]' < occupancy)"
if [[ "$occ" == "1.000" ]]; then
  echo "occupancy = 1"
else
  echo "occupancy not 1"
  exit 1
fi

cd ../../..





########################################
# Delete nt
########################################

number=${#nt_delete}
echo $number
if [ "$number" == 1 ]
then
  sed "/ $chain_delete   $nt_delete  /d" $pdb_file > PDB_without_nt/$pdb_id\_$chain_delete\_$nt_delete/omit.pdb
elif [ "$number" == 2 ]
then
  sed "/ $chain_delete  $nt_delete /d" $pdb_file > PDB_without_nt/$pdb_id\_$chain_delete\_$nt_delete/omit.pdb
elif [ "$number" == 3 ]
then
  sed "/ $chain_delete $nt_delete /d" $pdb_file > PDB_without_nt/$pdb_id\_$chain_delete\_$nt_delete/omit.pdb
elif [ "$number" == "4" ]
then
  sed "/ $chain_delete$nt_delete /d" $pdb_file > PDB_without_nt/$pdb_id\_$chain_delete\_$nt_delete/omit.pdb
fi


########################################
  # Run Phenix refine omit structure
  # Copy omit map to WC and HG folders
########################################
cd PDB_without_nt/$pdb_id\_$chain_delete\_$nt_delete
phenix.ready_set omit.pdb


# Check for omit.ligands.cif
FILE="omit.ligands.cif"
REFINE_SUCCESS=false


if [ -f "$FILE" ]; then
    echo "omit.updated.pdb $mtz_file omit.ligands.cif main.number_of_macro_cycles=5" >> phenix.refine.txt
    if phenix.refine omit.updated.pdb "$mtz_file" omit.ligands.cif \
        strategy=individual_sites+individual_adp+occupancies \
        main.number_of_macro_cycles=5; then
        
        if [ -f omit.updated_refine_001.mtz ]; then
            REFINE_SUCCESS=true
            cp omit.ligands.cif WC/
            cp omit.ligands.cif HG/
            cp omit.updated_refine_001.mtz WC/
            cp omit.updated_refine_001.mtz HG/
        fi
    fi
    
    if [ "$REFINE_SUCCESS" = false ]; then
        if phenix.refine omit.updated.pdb "$mtz_file" omit.ligands.cif \
            strategy=individual_sites+individual_adp+occupancies \
            main.number_of_macro_cycles=5 \
            xray_data.r_free_flags.generate=True overwrite=true; then
            
            if [ -f omit.updated_refine_001.mtz ]; then
                REFINE_SUCCESS=true
                echo "R-free generated" >> Rfactor_report.txt
                cp omit.ligands.cif WC/
                cp omit.ligands.cif HG/
                cp omit.updated_refine_001.mtz WC/
                cp omit.updated_refine_001.mtz HG/
            else
                echo "Refinement failed to produce MTZ file" >> Rfactor_report.txt
                exit 1
            fi
        else
            echo " Refinement failed" >> Rfactor_report.txt
            exit 1
        fi
    fi
    
else
    echo "omit.pdb $mtz_file main.number_of_macro_cycles=5" >> phenix.refine.txt
    
    if phenix.refine omit.pdb "$mtz_file" \
        strategy=individual_sites+individual_adp+occupancies \
        main.number_of_macro_cycles=5; then
        
        # Check if refinement completed successfully
        if [ -f omit_refine_001.mtz ]; then
            REFINE_SUCCESS=true

            echo "Refinement without nt" >> Rfactor_report.txt
            cp omit_refine_001.mtz WC/
            cp omit_refine_001.mtz HG/
        fi
    fi
    
    # If didn't complete, try with R-free flag generation
    if [ "$REFINE_SUCCESS" = false ]; then
        if phenix.refine omit.pdb "$mtz_file" \
            strategy=individual_sites+individual_adp+occupancies \
            main.number_of_macro_cycles=5 \
            xray_data.r_free_flags.generate=True overwrite=true; then
            
            if [ -f omit_refine_001.mtz ]; then
                REFINE_SUCCESS=true
                echo "R-free generated" >> Rfactor_report.txt
                cp omit_refine_001.mtz WC/
                cp omit_refine_001.mtz HG/
            else
                echo "Refinement failed to produce MTZ file" >> Rfactor_report.txt
                exit 1
            fi
        else
            echo " Refinement failed" >> Rfactor_report.txt
            exit 1
        fi
    fi
fi

########################################
  # Refinment in WC folder with omit map
########################################
cd WC
phenix.ready_set ${pdb_id}_final.pdb
FILE=${pdb_id}_final.ligands.cif
if [ -f "$FILE" ]
then
  phenix.refine ${pdb_id}_final.pdb omit.updated_refine_001.mtz ${pdb_id}_final.ligands.cif strategy=individual_sites+individual_adp+occupancies
  echo 'Done'
else
  phenix.refine ${pdb_id}_final.pdb omit_refine_001.mtz strategy=individual_sites+individual_adp+occupancies
  echo "Done 2"
fi

phenix.mtz2map ${pdb_id}_final_refine_001.mtz ${pdb_id}_final_refine_001.pdb


########################################
  # Make figure
########################################
script_positive="
from pymol import cmd
from pymol.cgo import COLOR, SPHERE
import math
import os
import sys
sys.path.append('/mnt/hdd_04/sg4213/Hoog-finder-2025/pymolprobity')
from pymolprobity import *


chain_1 = \"${chain_1}\"
nt_number_1 = \"${nt_number_1}\"
chain_2 = \"${chain_2}\"
nt_number_2 = \"${nt_number_2}\"
pdb_id = \"${pdb_id}\"

cmd.load(f\"{pdb_id}_final_refine_001.pdb\", \"pdb\")
cmd.load(f\"{pdb_id}_final_refine_001_2mFo-DFc.ccp4\", \"map\")
cmd.load(f\"{pdb_id}_final_refine_001_mFo-DFc.ccp4\", \"diffmap\")

#cmd.remove(\"solvent\")
cmd.select(\"solvent\")
cmd.extract(\"sol\", \"solvent\")
cmd.bg_color(\"white\")
cmd.color(\"lightblue\")

# Set cartoon ring mode (using a value, here 1)
cmd.set(\"cartoon_ring_mode\", 1)

# Select nucleic acids (for example DNA/RNA) and color them
cmd.select(\"DNA\", \"polymer.nucleic\")
cmd.color(\"gray50\", \"DNA\")

cmd.select(\"HG_pair\", f\"(chain {chain_1} and resi {nt_number_1}) or (chain {chain_2} and resi {nt_number_2})\")
cmd.extract(\"HG\", \"HG_pair\")
cmd.show(\"sticks\", \"HG\")

# Create isomeshes from the density maps around the 'HG' selection
cmd.isomesh(\"2mmFo-DFc\", \"map\", 2, \"HG\", carve=2)
cmd.isomesh(\"mmFo-DFc\", \"diffmap\", 3.0, \"HG\", carve=3)

# Color the density map objects for clarity
cmd.color(\"actinium\", \"2mmFo-DFc\")
cmd.color(\"barium\", \"mmFo-DFc\")

# Set mesh properties
cmd.set(\"mesh_negative_visible\", \"on\")
cmd.set(\"mesh_negative_color\", \"bohrium\")
cmd.set(\"mesh_negative_visible\", \"off\", \"2mmFo-DFc\")
cmd.set(\"mesh_width\", 0.4)

# Set rendering options
cmd.set(\"ray_trace_fog\", 0)
cmd.set(\"depth_cue\", 0)
cmd.set(\"ray_shadows\", \"off\")

# Color nitrogen and oxygen atoms
cmd.color(\"blue\", \"name N*\")
cmd.color(\"gold\", \"name O*\")

# Add hydrogens and create contact distance object
cmd.h_add(\"HG\")
main.reduce_object(\"HG\")
main.probe_object(\"HG\")

# Focus the view on specific atoms (here, those named C4', C5, and N3 in HG)
cmd.select(\"view\", \"name C4'+C5+N3 and HG\")
cmd.orient(\"view\")
cmd.zoom(\"view\")
cmd.center(\"view\")

# Remove the originally loaded pdb object to clear the workspace
cmd.delete(\"pdb\")

# Render the scene as a PNG image and then exit
cmd.png(f\"{pdb_id}_{chain_1}_{nt_number_1}_{chain_2}_{nt_number_2}_WC_map_water.png\")
cmd.turn(axis='x', angle=90.0)
cmd.set(\"ray_trace_frames\", 0)
cmd.png(f\"{pdb_id}_{chain_1}_{nt_number_1}_{chain_2}_{nt_number_2}_WC_map_water_90.png\")
cmd.quit() "


script_negative_both="

from pymol import cmd
from pymol.cgo import COLOR, SPHERE
import math
import os
import sys
sys.path.append('/mnt/hdd_04/sg4213/Hoog-finder-2025/pymolprobity')
from pymolprobity import *

chain_1 = \"${chain_1}\"
nt_number_1 = \"${nt_number_1}\"
chain_2 = \"${chain_2}\"
nt_number_2 = \"${nt_number_2}\"
pdb_id = \"${pdb_id}\"

cmd.load(f\"{pdb_id}_final_refine_001.pdb\", \"pdb\")
cmd.load(f\"{pdb_id}_final_refine_001_2mFo-DFc.ccp4\", \"map\")
cmd.load(f\"{pdb_id}_final_refine_001_mFo-DFc.ccp4\", \"diffmap\")

#cmd.remove(\"solvent\")
cmd.select(\"solvent\")
cmd.extract(\"sol\", \"solvent\")
cmd.bg_color(\"white\")
cmd.color(\"lightblue\")

# Set cartoon ring mode (using a value, here 1)
cmd.set(\"cartoon_ring_mode\", 1)

# Select nucleic acids (for example DNA/RNA) and color them
cmd.select(\"DNA\", \"polymer.nucleic\")
cmd.color(\"gray50\", \"DNA\")

cmd.select(\"HG_pair\", f\"(chain {chain_1} and resi \\{nt_number_1}) or (chain {chain_2} and resi \\{nt_number_2})\")
cmd.extract(\"HG\", \"HG_pair\")
cmd.show(\"sticks\", \"HG\")

# Create isomeshes from the density maps around the 'HG' selection
cmd.isomesh(\"2mmFo-DFc\", \"map\", 2, \"HG\", carve=2)
cmd.isomesh(\"mmFo-DFc\", \"diffmap\", 3.0, \"HG\", carve=3)

# Color the density map objects for clarity
cmd.color(\"actinium\", \"2mmFo-DFc\")
cmd.color(\"barium\", \"mmFo-DFc\")

# Set mesh properties
cmd.set(\"mesh_negative_visible\", \"on\")
cmd.set(\"mesh_negative_color\", \"bohrium\")
cmd.set(\"mesh_negative_visible\", \"off\", \"2mmFo-DFc\")
cmd.set(\"mesh_width\", 0.4)

# Set rendering options
cmd.set(\"ray_trace_fog\", 0)
cmd.set(\"depth_cue\", 0)
cmd.set(\"ray_shadows\", \"off\")

# Color nitrogen and oxygen atoms
cmd.color(\"blue\", \"name N*\")
cmd.color(\"gold\", \"name O*\")

# Add hydrogens and create contact distance object
cmd.h_add(\"HG\")
main.reduce_object(\"HG\")
main.probe_object(\"HG\")

# Focus the view on specific atoms (here, those named C4', C5, and N3 in HG)
cmd.select(\"view\", \"name C4'+C5+N3 and HG\")
cmd.orient(\"view\")
cmd.zoom(\"view\")
cmd.center(\"view\")

# Remove the originally loaded pdb object to clear the workspace
cmd.delete(\"pdb\")

# Render the scene as a PNG image and then exit
# Render the scene as a PNG image and then exit
cmd.png(f\"{pdb_id}_{chain_1}_{nt_number_1}_{chain_2}_{nt_number_2}_WC_map_water.png\")
cmd.turn(axis='x', angle=90.0)
cmd.set(\"ray_trace_frames\", 0)
cmd.png(f\"{pdb_id}_{chain_1}_{nt_number_1}_{chain_2}_{nt_number_2}_WC_map_water_90.png\")

cmd.quit() "

script_negative_1="

from pymol import cmd
from pymol.cgo import COLOR, SPHERE
import math
import os
import sys
sys.path.append('/mnt/hdd_04/sg4213/Hoog-finder-2025/pymolprobity')
from pymolprobity import *

chain_1 = \"${chain_1}\"
nt_number_1 = \"${nt_number_1}\"
chain_2 = \"${chain_2}\"
nt_number_2 = \"${nt_number_2}\"
pdb_id = \"${pdb_id}\"

cmd.load(f\"{pdb_id}_final_refine_001.pdb\", \"pdb\")
cmd.load(f\"{pdb_id}_final_refine_001_2mFo-DFc.ccp4\", \"map\")
cmd.load(f\"{pdb_id}_final_refine_001_mFo-DFc.ccp4\", \"diffmap\")

#cmd.remove(\"solvent\")
cmd.select(\"solvent\")   
cmd.extract(\"sol\", \"solvent\")
cmd.bg_color(\"white\")
cmd.color(\"lightblue\")

# Set cartoon ring mode (using a value, here 1)
cmd.set(\"cartoon_ring_mode\", 1)

# Select nucleic acids (for example DNA/RNA) and color them
cmd.select(\"DNA\", \"polymer.nucleic\")
cmd.color(\"gray50\", \"DNA\")

cmd.select(\"HG_pair\", f\"(chain {chain_1} and resi \\{nt_number_1}) or (chain {chain_2} and resi {nt_number_2})\")
cmd.extract(\"HG\", \"HG_pair\")
cmd.show(\"sticks\", \"HG\")

# Create isomeshes from the density maps around the 'HG' selection
cmd.isomesh(\"2mmFo-DFc\", \"map\", 2, \"HG\", carve=2)
cmd.isomesh(\"mmFo-DFc\", \"diffmap\", 3.0, \"HG\", carve=3)

# Color the density map objects for clarity
cmd.color(\"actinium\", \"2mmFo-DFc\")
cmd.color(\"barium\", \"mmFo-DFc\")

# Set mesh properties
cmd.set(\"mesh_negative_visible\", \"on\")
cmd.set(\"mesh_negative_color\", \"bohrium\")
cmd.set(\"mesh_negative_visible\", \"off\", \"2mmFo-DFc\")
cmd.set(\"mesh_width\", 0.4)

# Set rendering options
cmd.set(\"ray_trace_fog\", 0)
cmd.set(\"depth_cue\", 0)
cmd.set(\"ray_shadows\", \"off\")

# Color nitrogen and oxygen atoms
cmd.color(\"blue\", \"name N*\")
cmd.color(\"gold\", \"name O*\")

# Add hydrogens and create contact distance object
cmd.h_add(\"HG\")
main.reduce_object(\"HG\")
main.probe_object(\"HG\")

# Focus the view on specific atoms (here, those named C4', C5, and N3 in HG)
cmd.select(\"view\", \"name C4'+C5+N3 and HG\")
cmd.orient(\"view\")
cmd.zoom(\"view\")
cmd.center(\"view\")

# Remove the originally loaded pdb object to clear the workspace
cmd.delete(\"pdb\")

# Render the scene as a PNG image and then exit
# Render the scene as a PNG image and then exit
cmd.png(f\"{pdb_id}_{chain_1}_{nt_number_1}_{chain_2}_{nt_number_2}_WC_map_water.png\")
cmd.turn(axis='x', angle=90.0)
cmd.set(\"ray_trace_frames\", 0)
cmd.png(f\"{pdb_id}_{chain_1}_{nt_number_1}_{chain_2}_{nt_number_2}_WC_map_water_90.png\")
cmd.quit() "

script_negative_2="

from pymol import cmd
from pymol.cgo import COLOR, SPHERE
import math
import os
import sys
sys.path.append('/mnt/hdd_04/sg4213/Hoog-finder-2025/pymolprobity')
from pymolprobity import *

chain_1 = \"${chain_1}\"
nt_number_1 = \"${nt_number_1}\"
chain_2 = \"${chain_2}\"
nt_number_2 = \"${nt_number_2}\"
pdb_id = \"${pdb_id}\"

cmd.load(f\"{pdb_id}_final_refine_001.pdb\", \"pdb\")
cmd.load(f\"{pdb_id}_final_refine_001_2mFo-DFc.ccp4\", \"map\")
cmd.load(f\"{pdb_id}_final_refine_001_mFo-DFc.ccp4\", \"diffmap\")

#cmd.remove(\"solvent\")
cmd.select(\"solvent\")   
cmd.extract(\"sol\", \"solvent\")
cmd.bg_color(\"white\")
cmd.color(\"lightblue\")

# Set cartoon ring mode (using a value, here 1)
cmd.set(\"cartoon_ring_mode\", 1)

# Select nucleic acids (for example DNA/RNA) and color them
cmd.select(\"DNA\", \"polymer.nucleic\")
cmd.color(\"gray50\", \"DNA\")

cmd.select(\"HG_pair\", f\"(chain {chain_1} and resi {nt_number_1}) or (chain {chain_2} and resi \\{nt_number_2})\")
cmd.extract(\"HG\", \"HG_pair\")
cmd.show(\"sticks\", \"HG\")

# Create isomeshes from the density maps around the 'HG' selection
cmd.isomesh(\"2mmFo-DFc\", \"map\", 2, \"HG\", carve=2)
cmd.isomesh(\"mmFo-DFc\", \"diffmap\", 3.0, \"HG\", carve=3)

# Color the density map objects for clarity
cmd.color(\"actinium\", \"2mmFo-DFc\")
cmd.color(\"barium\", \"mmFo-DFc\")

# Set mesh properties
cmd.set(\"mesh_negative_visible\", \"on\")
cmd.set(\"mesh_negative_color\", \"bohrium\")
cmd.set(\"mesh_negative_visible\", \"off\", \"2mmFo-DFc\")
cmd.set(\"mesh_width\", 0.4)

# Set rendering options
cmd.set(\"ray_trace_fog\", 0)
cmd.set(\"depth_cue\", 0)
cmd.set(\"ray_shadows\", \"off\")

# Color nitrogen and oxygen atoms
cmd.color(\"blue\", \"name N*\")
cmd.color(\"gold\", \"name O*\")

# Add hydrogens and create contact distance object
cmd.h_add(\"HG\")
main.reduce_object(\"HG\")
main.probe_object(\"HG\")

# Focus the view on specific atoms (here, those named C4', C5, and N3 in HG)
cmd.select(\"view\", \"name C4'+C5+N3 and HG\")
cmd.orient(\"view\")
cmd.zoom(\"view\")
cmd.center(\"view\")

# Remove the originally loaded pdb object to clear the workspace
cmd.delete(\"pdb\")

# Render the scene as a PNG image and then exit
# Render the scene as a PNG image and then exit
cmd.png(f\"{pdb_id}_{chain_1}_{nt_number_1}_{chain_2}_{nt_number_2}_WC_map_water.png\")
cmd.turn(axis='x', angle=90.0)
cmd.set(\"ray_trace_frames\", 0)
cmd.png(f\"{pdb_id}_{chain_1}_{nt_number_1}_{chain_2}_{nt_number_2}_WC_map_water_90.png\")
cmd.quit() "


if [[ ${nt_number_1} -lt 0 && ${nt_number_2} -lt 0 ]]
then
    echo "$script_negative_both" > plot.py
elif [[ ${nt_number_1} -lt 0 && ${nt_number_2} -ge 0 ]]
then
    echo "$script_negative_1" > plot.py
elif [[ ${nt_number_2} -lt 0 && ${nt_number_1} -ge 0 ]]
then
    echo "$script_negative_2" > plot.py
elif [[ ${nt_number_2} -ge 0 && ${nt_number_1} -ge 0 ]]
then
    echo "$script_positive" > plot.py
else
    echo "$script_positive" > plot.py
fi
python3 plot.py
echo "PLOTTING"
cd ..
########################################
  # Refinment in HG folder with omit map
########################################
cd HG

echo flip.py ${pdb_id}_final $chain_delete $nt_delete
if [ $nt_delete -lt 0 ]
then
    python3 flip_negative.py ${pdb_id}_final $chain_delete $nt_delete
elif [ $nt_delete -ge 0 ]
then
    python3 flip.py ${pdb_id}_final $chain_delete $nt_delete
else
    python3 flip.py ${pdb_id}_final $chain_delete $nt_delete
fi



phenix.ready_set ${pdb_id}_final_flipped.pdb
FILE=${pdb_id}_final_flipped.ligands.cif
if [ -f "$FILE" ]
then
  echo ${pdb_id}_flipped.pdb omit.updated_refine_001.mtz ${pdb_id}_flipped.ligands.cif refine.sites.individual="(chain $chain_delete and resseq "$(($nt_delete-2))":"$(($nt_delete+2))")">> phenix.refine.txt
  phenix.refine ${pdb_id}_final_flipped.pdb omit.updated_refine_001.mtz ${pdb_id}_final_flipped.ligands.cif refine.sites.individual="(chain $chain_delete and resseq "$(($nt_delete-2))":"$(($nt_delete+2))" and not resname HOH)" miller_array.labels.name=F-obs

else
  echo ${pdb_id}_flipped.pdb *.mtz refine.sites.individual="(chain $chain_delete and resseq "$(($nt_delete-2))":"$(($nt_delete+2))")">> phenix.refine.txt
  phenix.refine ${pdb_id}_final_flipped.pdb *.mtz refine.sites.individual="(chain $chain_delete and resseq "$(($nt_delete-2))":"$(($nt_delete+2))" and not resname HOH)" miller_array.labels.name=F-obs
fi

FILE=${pdb_id}_final_flipped.ligands.cif
if [ -f "$FILE" ]
then
  echo ${pdb_id}_flipped_refine_001.pdb omit.updated_refine_001.mtz ${pdb_id}_flipped.ligands.cif >> phenix.refine.txt
  phenix.refine ${pdb_id}_final_flipped_refine_001.pdb omit.updated_refine_001.mtz ${pdb_id}_final_flipped.ligands.cif strategy=individual_sites+individual_adp+occupancies
  echo "Refinment in WC" >> Rfactor_report.txt

else
  echo ${pdb_id}_flipped_refine_001.pdb omit_refine_001.mtz >> phenix.refine.txt
  phenix.refine ${pdb_id}_final_flipped_refine_001.pdb omit_refine_001.mtz strategy=individual_sites+individual_adp+occupancies
  echo "Refinment without nt" >> Rfactor_report.txt
fi
phenix.mtz2map ${pdb_id}_final_flipped_refine_001_refine_001.mtz ${pdb_id}_final_flipped_refine_001_refine_001.pdb

########################################
  # Make figure
########################################
script_positive_hg="
from pymol import cmd
from pymol.cgo import COLOR, SPHERE
import math
import os
import sys
sys.path.append('/mnt/hdd_04/sg4213/Hoog-finder-2025/pymolprobity')
from pymolprobity import *

chain_1 = \"${chain_1}\"
nt_number_1 = \"${nt_number_1}\"
chain_2 = \"${chain_2}\"
nt_number_2 = \"${nt_number_2}\"
pdb_id = \"${pdb_id}\"

cmd.load(f\"{pdb_id}_final_flipped_refine_001_refine_001.pdb\", \"pdb\")
cmd.load(f\"{pdb_id}_final_flipped_refine_001_refine_001_2mFo-DFc.ccp4\", \"map\")
cmd.load(f\"{pdb_id}_final_flipped_refine_001_refine_001_mFo-DFc.ccp4\", \"diffmap\")

#cmd.remove(\"solvent\")
cmd.select(\"solvent\")   
cmd.extract(\"sol\", \"solvent\")
cmd.bg_color(\"white\")
cmd.color(\"lightblue\")

# Set cartoon ring mode (using a value, here 1)
cmd.set(\"cartoon_ring_mode\", 1)

# Select nucleic acids (for example DNA/RNA) and color them
cmd.select(\"DNA\", \"polymer.nucleic\")
cmd.color(\"gray50\", \"DNA\")

cmd.select(\"HG_pair\", f\"(chain {chain_1} and resi {nt_number_1}) or (chain {chain_2} and resi {nt_number_2})\")
cmd.extract(\"HG\", \"HG_pair\")
cmd.show(\"sticks\", \"HG\")

# Create isomeshes from the density maps around the 'HG' selection
cmd.isomesh(\"2mmFo-DFc\", \"map\", 2, \"HG\", carve=2)
cmd.isomesh(\"mmFo-DFc\", \"diffmap\", 3.0, \"HG\", carve=3)

# Color the density map objects for clarity
cmd.color(\"actinium\", \"2mmFo-DFc\")
cmd.color(\"barium\", \"mmFo-DFc\")

# Set mesh properties
cmd.set(\"mesh_negative_visible\", \"on\")
cmd.set(\"mesh_negative_color\", \"bohrium\")
cmd.set(\"mesh_negative_visible\", \"off\", \"2mmFo-DFc\")
cmd.set(\"mesh_width\", 0.4)

# Set rendering options
cmd.set(\"ray_trace_fog\", 0)
cmd.set(\"depth_cue\", 0)
cmd.set(\"ray_shadows\", \"off\")

# Color nitrogen and oxygen atoms
cmd.color(\"blue\", \"name N*\")
cmd.color(\"gold\", \"name O*\")

# Add hydrogens and create contact distance object
cmd.h_add(\"HG\")
main.reduce_object(\"HG\")
main.probe_object(\"HG\")

# Focus the view on specific atoms (here, those named C4', C5, and N3 in HG)
cmd.select(\"view\", \"name C4'+C5+N3 and HG\")
cmd.orient(\"view\")
cmd.zoom(\"view\")
cmd.center(\"view\")

# Remove the originally loaded pdb object to clear the workspace
cmd.delete(\"pdb\")

# Render the scene as a PNG image and then exit
cmd.png(f\"{pdb_id}_{chain_1}_{nt_number_1}_{chain_2}_{nt_number_2}_HG_map_water.png\")
cmd.turn(axis='x', angle=90.0)
cmd.set(\"ray_trace_frames\", 0)
cmd.png(f\"{pdb_id}_{chain_1}_{nt_number_1}_{chain_2}_{nt_number_2}_HG_map_water_90.png\")
cmd.quit() "


script_negative_both_hg="

from pymol import cmd
from pymol.cgo import COLOR, SPHERE
import math
import os
import sys
sys.path.append('/mnt/hdd_04/sg4213/Hoog-finder-2025/pymolprobity')
from pymolprobity import *


chain_1 = \"${chain_1}\"
nt_number_1 = \"${nt_number_1}\"
chain_2 = \"${chain_2}\"
nt_number_2 = \"${nt_number_2}\"
pdb_id = \"${pdb_id}\"

cmd.load(f\"{pdb_id}_final_flipped_refine_001_refine_001.pdb\", \"pdb\")
cmd.load(f\"{pdb_id}_final_flipped_refine_001_refine_001_2mFo-DFc.ccp4\", \"map\")
cmd.load(f\"{pdb_id}_final_flipped_refine_001_refine_001_mFo-DFc.ccp4\", \"diffmap\")

#cmd.remove(\"solvent\")
cmd.select(\"solvent\")   
cmd.extract(\"sol\", \"solvent\")
cmd.bg_color(\"white\")
cmd.color(\"lightblue\")

# Set cartoon ring mode (using a value, here 1)
cmd.set(\"cartoon_ring_mode\", 1)

# Select nucleic acids (for example DNA/RNA) and color them
cmd.select(\"DNA\", \"polymer.nucleic\")
cmd.color(\"gray50\", \"DNA\")

cmd.select(\"HG_pair\", f\"(chain {chain_1} and resi \\{nt_number_1}) or (chain {chain_2} and resi \\{nt_number_2})\")
cmd.extract(\"HG\", \"HG_pair\")
cmd.show(\"sticks\", \"HG\")

# Create isomeshes from the density maps around the 'HG' selection
cmd.isomesh(\"2mmFo-DFc\", \"map\", 2, \"HG\", carve=2)
cmd.isomesh(\"mmFo-DFc\", \"diffmap\", 3.0, \"HG\", carve=3)

# Color the density map objects for clarity
cmd.color(\"actinium\", \"2mmFo-DFc\")
cmd.color(\"barium\", \"mmFo-DFc\")

# Set mesh properties
cmd.set(\"mesh_negative_visible\", \"on\")
cmd.set(\"mesh_negative_color\", \"bohrium\")
cmd.set(\"mesh_negative_visible\", \"off\", \"2mmFo-DFc\")
cmd.set(\"mesh_width\", 0.4)

# Set rendering options
cmd.set(\"ray_trace_fog\", 0)
cmd.set(\"depth_cue\", 0)
cmd.set(\"ray_shadows\", \"off\")

# Color nitrogen and oxygen atoms
cmd.color(\"blue\", \"name N*\")
cmd.color(\"gold\", \"name O*\")

# Add hydrogens and create contact distance object
cmd.h_add(\"HG\")
main.reduce_object(\"HG\")
main.probe_object(\"HG\")

# Focus the view on specific atoms (here, those named C4', C5, and N3 in HG)
cmd.select(\"view\", \"name C4'+C5+N3 and HG\")
cmd.orient(\"view\")
cmd.zoom(\"view\")
cmd.center(\"view\")

# Remove the originally loaded pdb object to clear the workspace
cmd.delete(\"pdb\")

# Render the scene as a PNG image and then exit
cmd.png(f\"{pdb_id}_{chain_1}_{nt_number_1}_{chain_2}_{nt_number_2}_HG_map_water.png\")
cmd.turn(axis='x', angle=90.0)
cmd.set(\"ray_trace_frames\", 0)
cmd.png(f\"{pdb_id}_{chain_1}_{nt_number_1}_{chain_2}_{nt_number_2}_HG_map_water_90.png\")
cmd.quit() "

script_negative_1_hg="

from pymol import cmd
from pymol.cgo import COLOR, SPHERE
import math
import os
import sys
sys.path.append('/mnt/hdd_04/sg4213/Hoog-finder-2025/pymolprobity')
from pymolprobity import *


chain_1 = \"${chain_1}\"
nt_number_1 = \"${nt_number_1}\"
chain_2 = \"${chain_2}\"
nt_number_2 = \"${nt_number_2}\"
pdb_id = \"${pdb_id}\"

cmd.load(f\"{pdb_id}_final_flipped_refine_001_refine_001.pdb\", \"pdb\")
cmd.load(f\"{pdb_id}_final_flipped_refine_001_refine_001_2mFo-DFc.ccp4\", \"map\")
cmd.load(f\"{pdb_id}_final_flipped_refine_001_refine_001_mFo-DFc.ccp4\", \"diffmap\")

#cmd.remove(\"solvent\")
cmd.select(\"solvent\")   
cmd.extract(\"sol\", \"solvent\")
cmd.bg_color(\"white\")
cmd.color(\"lightblue\")

# Set cartoon ring mode (using a value, here 1)
cmd.set(\"cartoon_ring_mode\", 1)

# Select nucleic acids (for example DNA/RNA) and color them
cmd.select(\"DNA\", \"polymer.nucleic\")
cmd.color(\"gray50\", \"DNA\")

cmd.select(\"HG_pair\", f\"(chain {chain_1} and resi \\{nt_number_1}) or (chain {chain_2} and resi {nt_number_2})\")
cmd.extract(\"HG\", \"HG_pair\")
cmd.show(\"sticks\", \"HG\")

# Create isomeshes from the density maps around the 'HG' selection
cmd.isomesh(\"2mmFo-DFc\", \"map\", 2, \"HG\", carve=2)
cmd.isomesh(\"mmFo-DFc\", \"diffmap\", 3.0, \"HG\", carve=3)

# Color the density map objects for clarity
cmd.color(\"actinium\", \"2mmFo-DFc\")
cmd.color(\"barium\", \"mmFo-DFc\")

# Set mesh properties
cmd.set(\"mesh_negative_visible\", \"on\")
cmd.set(\"mesh_negative_color\", \"bohrium\")
cmd.set(\"mesh_negative_visible\", \"off\", \"2mmFo-DFc\")
cmd.set(\"mesh_width\", 0.4)

# Set rendering options
cmd.set(\"ray_trace_fog\", 0)
cmd.set(\"depth_cue\", 0)
cmd.set(\"ray_shadows\", \"off\")

# Color nitrogen and oxygen atoms
cmd.color(\"blue\", \"name N*\")
cmd.color(\"gold\", \"name O*\")

# Add hydrogens and create contact distance object
cmd.h_add(\"HG\")
main.reduce_object(\"HG\")
main.probe_object(\"HG\")

# Focus the view on specific atoms (here, those named C4', C5, and N3 in HG)
cmd.select(\"view\", \"name C4'+C5+N3 and HG\")
cmd.orient(\"view\")
cmd.zoom(\"view\")
cmd.center(\"view\")

# Remove the originally loaded pdb object to clear the workspace
cmd.delete(\"pdb\")

# Render the scene as a PNG image and then exit
cmd.png(f\"{pdb_id}_{chain_1}_{nt_number_1}_{chain_2}_{nt_number_2}_HG_map_water.png\")
cmd.turn(axis='x', angle=90.0)
cmd.set(\"ray_trace_frames\", 0)
cmd.png(f\"{pdb_id}_{chain_1}_{nt_number_1}_{chain_2}_{nt_number_2}_HG_map_water_90.png\")
cmd.quit() "

script_negative_2_hg="

from pymol import cmd
from pymol.cgo import COLOR, SPHERE
import math
import os
import sys
sys.path.append('/mnt/hdd_04/sg4213/Hoog-finder-2025/pymolprobity')
from pymolprobity import *


chain_1 = \"${chain_1}\"
nt_number_1 = \"${nt_number_1}\"
chain_2 = \"${chain_2}\"
nt_number_2 = \"${nt_number_2}\"
pdb_id = \"${pdb_id}\"

cmd.load(f\"{pdb_id}_final_flipped_refine_001_refine_001.pdb\", \"pdb\")
cmd.load(f\"{pdb_id}_final_flipped_refine_001_refine_001_2mFo-DFc.ccp4\", \"map\")
cmd.load(f\"{pdb_id}_final_flipped_refine_001_refine_001_mFo-DFc.ccp4\", \"diffmap\")

#cmd.remove(\"solvent\")
cmd.select(\"solvent\")   
cmd.extract(\"sol\", \"solvent\")
cmd.bg_color(\"white\")
cmd.color(\"lightblue\")

# Set cartoon ring mode (using a value, here 1)
cmd.set(\"cartoon_ring_mode\", 1)

# Select nucleic acids (for example DNA/RNA) and color them
cmd.select(\"DNA\", \"polymer.nucleic\")
cmd.color(\"gray50\", \"DNA\")

cmd.select(\"HG_pair\", f\"(chain {chain_1} and resi {nt_number_1}) or (chain {chain_2} and resi \\{nt_number_2})\")
cmd.extract(\"HG\", \"HG_pair\")
cmd.show(\"sticks\", \"HG\")

# Create isomeshes from the density maps around the 'HG' selection
cmd.isomesh(\"2mmFo-DFc\", \"map\", 2, \"HG\", carve=2)
cmd.isomesh(\"mmFo-DFc\", \"diffmap\", 3.0, \"HG\", carve=3)

# Color the density map objects for clarity
cmd.color(\"actinium\", \"2mmFo-DFc\")
cmd.color(\"barium\", \"mmFo-DFc\")

# Set mesh properties
cmd.set(\"mesh_negative_visible\", \"on\")
cmd.set(\"mesh_negative_color\", \"bohrium\")
cmd.set(\"mesh_negative_visible\", \"off\", \"2mmFo-DFc\")
cmd.set(\"mesh_width\", 0.4)

# Set rendering options
cmd.set(\"ray_trace_fog\", 0)
cmd.set(\"depth_cue\", 0)
cmd.set(\"ray_shadows\", \"off\")

# Color nitrogen and oxygen atoms
cmd.color(\"blue\", \"name N*\")
cmd.color(\"gold\", \"name O*\")

# Add hydrogens and create contact distance object
cmd.h_add(\"HG\")
main.reduce_object(\"HG\")
main.probe_object(\"HG\")

# Focus the view on specific atoms (here, those named C4', C5, and N3 in HG)
cmd.select(\"view\", \"name C4'+C5+N3 and HG\")
cmd.orient(\"view\")
cmd.zoom(\"view\")
cmd.center(\"view\")

# Remove the originally loaded pdb object to clear the workspace
cmd.delete(\"pdb\")

# Render the scene as a PNG image and then exit
cmd.png(f\"{pdb_id}_{chain_1}_{nt_number_1}_{chain_2}_{nt_number_2}_HG_map_water.png\")
cmd.turn(axis='x', angle=90.0)
cmd.set(\"ray_trace_frames\", 0)
cmd.png(f\"{pdb_id}_{chain_1}_{nt_number_1}_{chain_2}_{nt_number_2}_HG_map_water_90.png\")
cmd.quit() "

if [[ ${nt_number_1} -lt 0 && ${nt_number_2} -lt 0 ]]
then
    echo "$script_negative_both_hg" > plot.py
elif [[ ${nt_number_1} -lt 0 && ${nt_number_2} -ge 0 ]]
then
    echo "$script_negative_1_hg" > plot.py
elif [[ ${nt_number_2} -lt 0 && ${nt_number_1} -ge 0 ]]
then
    echo "$script_negative_2_hg" > plot.py
elif [[ ${nt_number_2} -ge 0 && ${nt_number_1} -ge 0 ]]
then
    echo "$script_positive_hg" > plot.py
else
    echo "$script_positive_hg" > plot.py
fi

python3 plot.py
echo "PLOTTING"
cd ../..
