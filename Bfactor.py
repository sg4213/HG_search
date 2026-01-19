#!/usr/bin/env python3
import os
import sys
import glob
import numpy as np

# -------------------------
# Config
# -------------------------
ROOT_GLOBS = [
    "PDB_without_nt",
    "../PDB_without_nt",
]
HG_PDB_NAME = "{pdbid}_final_flipped_refine_001_refine_001.pdb"
WC_PDB_NAME = "{pdbid}_final_refine_001.pdb"

OUT_SUMMARY = "classification_files/Bfactor_summary.txt"

# -------------------------
# Helpers
# -------------------------
def get_purine_info(chain1, nt_type1, nt1, chain2, nt_type2, nt2):
    """
    Decide which side is purine; return (chain_purine, nt_purine).
    """
    if nt_type1 in ["A", "G"]:
        return chain1, nt1
    if nt_type2 in ["A", "G"]:
        return chain2, nt2
    return None, None

def locate_tuple_dir(pdbid, chain_purine, nt_purine):
    """
    Return the first tuple directory matching any ROOT_GLOBS:
      <root>/{pdbid}_{chain_purine}_{nt_purine}
    """
    for root in ROOT_GLOBS:
        tup = os.path.join(root, f"{pdbid}_{chain_purine}_{nt_purine}")
        if os.path.isdir(tup):
            return tup
    return None

def mean_b_for_residue(pdb_file: str, chain_id: str, resi: int) -> float | None:
    """
    Calculate mean B-factor for all atoms in a specific residue.
    """
    if not pdb_file or not os.path.exists(pdb_file):
        return None
    
    b_factors = []
    with open(pdb_file) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if line[21] != chain_id:
                continue
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            if resseq != resi:
                continue
            try:
                b = float(line[60:66])
                b_factors.append(b)
            except ValueError:
                pass
    
    if b_factors:
        mean_b = float(np.mean(b_factors))
        return round(mean_b, 1)  # precision: 1 digit after "."
    return None

# -------------------------
# Main
# -------------------------
if __name__ == "__main__":

    pdbid = sys.argv[1].lower()
    chain1 = sys.argv[2]
    nt_type1 = sys.argv[3]
    nt1 = int(sys.argv[4])
    chain2 = sys.argv[5]
    nt_type2 = sys.argv[6]
    nt2 = int(sys.argv[7])

    # Get purine information
    chain_purine, nt_purine = get_purine_info(
        chain1, nt_type1, nt1, chain2, nt_type2, nt2
    )
    
    if chain_purine is None:
        print(f" No purine found for {pdbid} {chain1}:{nt_type1}{nt1} - {chain2}:{nt_type2}{nt2}", file=sys.stderr)
        sys.exit(1)

    # Locate tuple directory
    tup_dir = locate_tuple_dir(pdbid, chain_purine, nt_purine)
    if not tup_dir:
        print(f" Missing tuple dir for {pdbid}_{chain_purine}_{nt_purine}", file=sys.stderr)
        # Write None values and exit
        os.makedirs("classification_files", exist_ok=True)
        if not os.path.exists(OUT_SUMMARY):
            with open(OUT_SUMMARY, "w") as f:
                f.write("pdb_id chain_1 nt_type_1 nt_number_1 chain_2 nt_type_2 nt_number_2 "
                        "mean_B_HG mean_B_WC\n")
        with open(OUT_SUMMARY, "a") as f:
            f.write(f"{pdbid} {chain1} {nt_type1} {nt1} {chain2} {nt_type2} {nt2} "
                    f"None None\n")
        sys.exit(1)


    # Locate PDB files
    hg_pdb = os.path.join(tup_dir, "HG", HG_PDB_NAME.format(pdbid=pdbid))
    wc_pdb = os.path.join(tup_dir, "WC", WC_PDB_NAME.format(pdbid=pdbid))

    # Calculate mean B-factors
    mean_b_hg = mean_b_for_residue(hg_pdb, chain_purine, nt_purine) if os.path.exists(hg_pdb) else None
    mean_b_wc = mean_b_for_residue(wc_pdb, chain_purine, nt_purine) if os.path.exists(wc_pdb) else None


    if not os.path.exists(OUT_SUMMARY):
        with open(OUT_SUMMARY, "w") as f:
            f.write("pdb_id chain_1 nt_type_1 nt_number_1 chain_2 nt_type_2 nt_number_2 "
                    "mean_B_HG mean_B_WC\n")

    # Write to summary file
    with open(OUT_SUMMARY, "a") as f:
        f.write(f"{pdbid} {chain1} {nt_type1} {nt1} {chain2} {nt_type2} {nt2} "
                f"{mean_b_hg} {mean_b_wc}\n")
