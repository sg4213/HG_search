#!/usr/bin/env python3
import os
import sys
import re
from pymol import cmd

# -------------------------
# Config
# -------------------------
ROOT_GLOBS = [
    "PDB_without_nt",
    "../PDB_without_nt",
]
HG_PDB_NAME = "{pdbid}_final_flipped_refine_001_refine_001.pdb"
WC_PDB_NAME = "{pdbid}_final_refine_001.pdb"

NEIGH_PDB_SUFFIX = "_bp_plus1.pdb"
NEIGH_TXT_SUFFIX = "_clashscore_local_neighbours.txt"
BP_PDB_SUFFIX = "_bp.pdb"
BP_TXT_SUFFIX = "_clashscore_local_bp.txt"
GLOBAL_TXT_SUFFIX = "_clashscore_global.txt"

OUT_SUMMARY = "classification_files/clashscore_summary.txt"

# -------------------------
# Helpers
# -------------------------
def get_purine_info(chain1, nt_type1, nt1, chain2, nt_type2, nt2):
    """
    Decide which side is purine; return (chain_purine, nt_purine, chain_pyrimidine, nt_pyrimidine).
    """
    if nt_type1 in ["A", "G"]:
        return chain1, nt1, chain2, nt2
    if nt_type2 in ["A", "G"]:
        return chain2, nt2, chain1, nt1
    return None, None, None, None

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

def extract_clashscore(txt_path):
    """Extract numeric clashscore from phenix.clashscore output file."""
    try:
        with open(txt_path, "r") as f:
            for line in f:
                if "clashscore" in line:
                    m = re.search(r"clashscore\s*[:=]\s*([0-9.]+)", line, flags=re.I)
                    if m:
                        return float(m.group(1))
    except Exception:
        pass
    return None

def run_clashscore(pdb_path, txt_out, force=False):
    """
    Run phenix.clashscore on pdb_path and write txt_out.
    If txt_out exists and force=False, skip recomputation.
    """
    if not os.path.exists(pdb_path):
        return False
    if os.path.exists(txt_out) and not force:
        return True
    rc = os.system(f"phenix.clashscore {pdb_path} > {txt_out} 2>&1")
    return rc == 0

def clamp_resi(n):
    return max(1, int(n))

def sel_neighbours(chain, nt):
    n = clamp_resi(nt)
    return (
        f"(chain {chain} and resi {n-1}) or "
        f"(chain {chain} and resi {n}) or "
        f"(chain {chain} and resi {n+1})"
    )

def sel_bp(chain1, nt1, chain2, nt2):
    return (
        f"(chain {chain1} and resi {int(nt1)}) or "
        f"(chain {chain2} and resi {int(nt2)})"
    )

def ensure_locals_for_state(tup_dir, state, pdb_name, name_id, chain_purine, nt_purine, chain_pyrimidine, nt_pyrimidine):
    """
    For given state dir (HG/WC):
      - load PDB
      - save two PDBs (neighbours & bp)
      - run clashscore on each → txts
    Return (neigh_txt_path, bp_txt_path), each may be None on failure.
    """
    state_dir = os.path.join(tup_dir, state)
    if not os.path.isdir(state_dir):
        print(f" Missing state dir: {state_dir}")
        return None, None

    pdb_path = os.path.join(state_dir, pdb_name)
    if not os.path.exists(pdb_path):
        print(f" Missing PDB: {pdb_path}")
        return None, None

    try:
        cmd.reinitialize()
        cmd.load(pdb_path, name_id)

        # Neighbours (nt_purine±1)
        cmd.select("sel_neigh", sel_neighbours(chain_purine, nt_purine))
        out_neigh_pdb = os.path.join(state_dir, f"{name_id}{NEIGH_PDB_SUFFIX}")
        cmd.save(out_neigh_pdb, "sel_neigh")
        out_neigh_txt = os.path.join(state_dir, f"{name_id}{NEIGH_TXT_SUFFIX}")
        run_clashscore(out_neigh_pdb, out_neigh_txt)

        # Base pair only
        cmd.select("sel_bp", sel_bp(chain_purine, nt_purine, chain_pyrimidine, nt_pyrimidine))
        out_bp_pdb = os.path.join(state_dir, f"{name_id}{BP_PDB_SUFFIX}")
        cmd.save(out_bp_pdb, "sel_bp")
        out_bp_txt = os.path.join(state_dir, f"{name_id}{BP_TXT_SUFFIX}")
        run_clashscore(out_bp_pdb, out_bp_txt)

        return out_neigh_txt, out_bp_txt

    except Exception as e:
        print(f" PyMOL error [{state}] {name_id}: {e}")
        return None, None

def ensure_global_for_state(tup_dir, state, pdb_name, name_id):
    """
    For given state dir (HG/WC):
      - run clashscore on the full model PDB
      - write *_clashscore_global.txt
    Return txt path or None on failure.
    """
    state_dir = os.path.join(tup_dir, state)
    if not os.path.isdir(state_dir):
        print(f" Missing state dir: {state_dir}")
        return None
    pdb_path = os.path.join(state_dir, pdb_name)
    if not os.path.exists(pdb_path):
        print(f" Missing PDB: {pdb_path}")
        return None
    out_txt = os.path.join(state_dir, f"{name_id}{GLOBAL_TXT_SUFFIX}")
    run_clashscore(pdb_path, out_txt)
    return out_txt if os.path.exists(out_txt) else None

# -------------------------
# Main
# -------------------------
if __name__ == "__main__":
    if len(sys.argv) != 8:
        print("Usage: python3 script.py pdb_id chain_1 nt_type_1 nt_number_1 chain_2 nt_type_2 nt_number_2", file=sys.stderr)
        sys.exit(1)

    pdbid = sys.argv[1]
    chain1 = sys.argv[2]
    nt_type1 = sys.argv[3]
    nt1 = int(sys.argv[4])
    chain2 = sys.argv[5]
    nt_type2 = sys.argv[6]
    nt2 = int(sys.argv[7])

    # Get purine information
    chain_purine, nt_purine, chain_pyrimidine, nt_pyrimidine = get_purine_info(
        chain1, nt_type1, nt1, chain2, nt_type2, nt2
    )
    
    if chain_purine is None:
        print(f" No purine found for {pdbid} {chain1}:{nt_type1}{nt1} - {chain2}:{nt_type2}{nt2}", file=sys.stderr)
        sys.exit(1)

    name_id = f"{pdbid}_{chain_purine}_{nt_purine}_{chain_pyrimidine}_{nt_pyrimidine}"

    tup_dir = locate_tuple_dir(pdbid, chain_purine, nt_purine)
    if not tup_dir:
        print(f" Missing tuple dir for {pdbid}_{chain_purine}_{nt_purine}", file=sys.stderr)
        # Write None values and exit
        os.makedirs("classification_files", exist_ok=True)
        if not os.path.exists(OUT_SUMMARY):
            with open(OUT_SUMMARY, "w") as f:
                f.write("pdb_id chain_1 nt_type_1 nt_number_1 chain_2 nt_type_2 nt_number_2 "
                        "WC_clashscore_global HG_clashscore_global "
                        "WC_clashscore_bp HG_clashscore_bp "
                        "WC_clashscore_neighbour HG_clashscore_neighbour\n")
        with open(OUT_SUMMARY, "a") as f:
            f.write(f"{pdbid} {chain1} {nt_type1} {nt1} {chain2} {nt_type2} {nt2} "
                    f"None None None None None None\n")
        sys.exit(1)


    # Process HG
    ensure_locals_for_state(
        tup_dir, "HG",
        HG_PDB_NAME.format(pdbid=pdbid),
        name_id, chain_purine, nt_purine, chain_pyrimidine, nt_pyrimidine,
    )
    ensure_global_for_state(
        tup_dir, "HG", HG_PDB_NAME.format(pdbid=pdbid), name_id
    )

    # Process WC
    ensure_locals_for_state(
        tup_dir, "WC",
        WC_PDB_NAME.format(pdbid=pdbid),
        name_id, chain_purine, nt_purine, chain_pyrimidine, nt_pyrimidine,
    )
    ensure_global_for_state(
        tup_dir, "WC", WC_PDB_NAME.format(pdbid=pdbid), name_id
    )

    # Collect scores
    hg_global_txt = os.path.join(tup_dir, "HG", f"{name_id}{GLOBAL_TXT_SUFFIX}")
    wc_global_txt = os.path.join(tup_dir, "WC", f"{name_id}{GLOBAL_TXT_SUFFIX}")
    hg_bp_txt = os.path.join(tup_dir, "HG", f"{name_id}{BP_TXT_SUFFIX}")
    wc_bp_txt = os.path.join(tup_dir, "WC", f"{name_id}{BP_TXT_SUFFIX}")
    hg_neigh_txt = os.path.join(tup_dir, "HG", f"{name_id}{NEIGH_TXT_SUFFIX}")
    wc_neigh_txt = os.path.join(tup_dir, "WC", f"{name_id}{NEIGH_TXT_SUFFIX}")

    wc_global = extract_clashscore(wc_global_txt) if os.path.exists(wc_global_txt) else None
    hg_global = extract_clashscore(hg_global_txt) if os.path.exists(hg_global_txt) else None
    wc_bp = extract_clashscore(wc_bp_txt) if os.path.exists(wc_bp_txt) else None
    hg_bp = extract_clashscore(hg_bp_txt) if os.path.exists(hg_bp_txt) else None
    wc_neigh = extract_clashscore(wc_neigh_txt) if os.path.exists(wc_neigh_txt) else None
    hg_neigh = extract_clashscore(hg_neigh_txt) if os.path.exists(hg_neigh_txt) else None

    # Create output directory and write header if needed
    os.makedirs("classification_files", exist_ok=True)
    if not os.path.exists(OUT_SUMMARY):
        with open(OUT_SUMMARY, "w") as f:
            f.write("pdb_id chain_1 nt_type_1 nt_number_1 chain_2 nt_type_2 nt_number_2 "
                    "WC_clashscore_global HG_clashscore_global "
                    "WC_clashscore_bp HG_clashscore_bp "
                    "WC_clashscore_neighbour HG_clashscore_neighbour\n")

    # Write to summary file
    with open(OUT_SUMMARY, "a") as f:
        f.write(f"{pdbid} {chain1} {nt_type1} {nt1} {chain2} {nt_type2} {nt2} "
                f"{wc_global} {hg_global} {wc_bp} {hg_bp} {wc_neigh} {hg_neigh}\n")