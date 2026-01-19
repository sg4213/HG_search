#!/usr/bin/env python3
import sys
import os
import math
import shutil
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# -------------------------
# Config
# -------------------------
COMBINED_TABLE = "classification_files/combined_metrics.txt"
PATH_TO_IMAGES = "PDB_without_nt"
OUTPUT_FOLDER = "reports"

# -------------------------
# Helper Functions
# -------------------------
def _is_none_or_nan(x) -> bool:
    """Check if value is None or NaN."""
    if x is None:
        return True
    if isinstance(x, str) and x.lower() == 'none':
        return True
    if isinstance(x, float) and math.isnan(x):
        return True
    return False

def _is_num(x):
    """Check if value is a valid number."""
    if _is_none_or_nan(x):
        return False
    try:
        float(x)
        return True
    except:
        return False

def interpret_metric(
    diff,
    metric_type: str = "default",
    tol_rscc: float = 0.007,
    tol_edia: float = 0.010,
    tol_other: float = 0.005,
):

    if _is_none_or_nan(diff):
        return None

    mt = metric_type.strip().lower().replace("_", "-")

    if mt == "rscc":
        tol = tol_rscc
        if -tol <= diff <= tol:
            return "Ambiguous"
        return "WC" if diff > tol else "HG"

    elif mt == "edia":
        tol = tol_edia
        if -tol <= diff <= tol:
            return "Ambiguous"
        return "WC" if diff > tol else "HG"

    elif mt == "clashscore":
        if diff == 0:
            return "Ambiguous"
        # lower clashscore is better
        return "HG" if diff > 0 else "WC"

    elif mt in ("r-work", "r-free", "rwork", "rfree"):
        return "Highlight in Red" if abs(diff) > 0.01 else "Normal"

    else:
        # default
        tol = tol_other
        if -tol <= diff <= tol:
            return "Ambiguous"
        return "WC" if diff > tol else "HG"

def determine_overall(result_dict):
    """
    Decide overall: return 'WC' if WC votes > HG votes, 'HG' if HG votes > WC votes,
    otherwise 'Ambiguous' (i.e., only when counts are equal).
    """
    wc_count = sum(v == "WC" for v in result_dict.values() if v is not None)
    hg_count = sum(v == "HG" for v in result_dict.values() if v is not None)
    
    if wc_count > hg_count:
        return "WC"
    if hg_count > wc_count:
        return "HG"
    return "Ambiguous"

def safe_float(val, default='NA'):
    """Convert value to float, return default if not possible."""
    if _is_none_or_nan(val):
        return default
    try:
        return float(val)
    except:
        return default

def format_value(val, decimals=3):
    """Format numeric value with specified decimals, or return 'NA'."""
    if _is_none_or_nan(val):
        return 'NA'
    try:
        f = float(val)
        return f"{f:.{decimals}f}"
    except:
        return 'NA'

def delete_folder_safely(folder_path):
    """Delete a folder and all its contents safely."""
    if os.path.exists(folder_path) and os.path.isdir(folder_path):
        try:
            shutil.rmtree(folder_path)
            return True
        except Exception as e:
            print(f"  Warning: Could not delete folder {folder_path}: {e}")
            return False
    else:
        print(f"  Folder does not exist: {folder_path}")
        return False

# -------------------------
# Main
# -------------------------
if __name__ == "__main__":
    if len(sys.argv) != 11:
        print("Usage: python3 make_report.py <pdb_code> <resolution> <chain_1> <nt_type_1> <nt_number_1> <chain_2> <nt_type_2> <nt_number_2> <chi_1> <chi_2>", file=sys.stderr)
        sys.exit(1)

    # Parse arguments
    pdb_code = sys.argv[1].lower()
    resolution = sys.argv[2]
    chain_1 = sys.argv[3]
    nt_type_1 = sys.argv[4].upper()
    nt_number_1 = int(sys.argv[5])
    chain_2 = sys.argv[6]
    nt_type_2 = sys.argv[7].upper()
    nt_number_2 = int(sys.argv[8])
    chi_1 = sys.argv[9]
    chi_2 = sys.argv[10]

    # Check if combined table exists
    if not os.path.exists(COMBINED_TABLE):
        print(f" Error: Combined metrics table not found: {COMBINED_TABLE}", file=sys.stderr)
        sys.exit(1)

    # Load the data
    try:
        data = pd.read_csv(COMBINED_TABLE, sep=r'\s+')
    except Exception as e:
        print(f" Error reading combined table: {e}", file=sys.stderr)
        sys.exit(1)

    # Find the specific entry
    mask = (
        (data['pdb_id'].str.lower() == pdb_code) &
        (data['chain_1'] == chain_1) &
        (data['nt_type_1'] == nt_type_1) &
        (data['nt_number_1'] == nt_number_1) &
        (data['chain_2'] == chain_2) &
        (data['nt_type_2'] == nt_type_2) &
        (data['nt_number_2'] == nt_number_2)
    )
    
    matching_rows = data[mask]
    
    if len(matching_rows) == 0:
        print(f" Error: No entry found for {pdb_code} {chain_1} {nt_type_1} {nt_number_1} {chain_2} {nt_type_2} {nt_number_2}")
        sys.exit(1)
    
    if len(matching_rows) > 1:
        print(f"  Warning: Multiple entries found, using the first one")
    
    row_index = matching_rows.index[0]
    row = matching_rows.iloc[0]

    # Create output directory
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    
    # Output PDF name
    OUTPUT_PDF = f"{OUTPUT_FOLDER}/{pdb_code}_{chain_1}_{nt_type_1}_{nt_number_1}_{chain_2}_{nt_type_2}_{nt_number_2}.pdf"

    try:
        # Determine purine and its conformation
        if nt_type_1 in ['A', 'G']:
            chain_purine = chain_1
            nt_purine = nt_number_1
            purine = nt_type_1
            conformation = chi_1
        elif nt_type_2 in ['A', 'G']:
            chain_purine = chain_2
            nt_purine = nt_number_2
            purine = nt_type_2
            conformation = chi_2
        else:
            print(f" Error: No purine found in base pair")
            sys.exit(1)

        # Determine base pair type from conformation
        if conformation == 'syn':
            conformation_bp = 'HG'
        elif conformation == 'anti':
            conformation_bp = 'WC'
        else:
            conformation_bp = 'Other'

        # Swap logic: if conformation is syn, swap WC/HG prefixes and image paths
        if conformation == 'syn':
            # Swap: columns labeled WC are the "HG" values to compare against,
            # and columns labeled HG are the "WC" values. Also swap image files under titles.
            hg_prefix, wc_prefix = 'WC', 'HG'
            hg_img_suffix, wc_img_suffix = 'WC_map_water', 'HG_map_water'
        else:
            # Default: use columns as labeled
            hg_prefix, wc_prefix = 'HG', 'WC'
            hg_img_suffix, wc_img_suffix = 'HG_map_water', 'WC_map_water'

        # Generate paths using swapped prefixes
        base_path = os.path.join(PATH_TO_IMAGES, f'{pdb_code}_{chain_purine}_{nt_purine}')
        hg_image = os.path.join(base_path, hg_prefix, f'{pdb_code}_{chain_1}_{nt_number_1}_{chain_2}_{nt_number_2}_{hg_img_suffix}.png')
        wc_image = os.path.join(base_path, wc_prefix, f'{pdb_code}_{chain_1}_{nt_number_1}_{chain_2}_{nt_number_2}_{wc_img_suffix}.png')
        hg_image_90 = os.path.join(base_path, hg_prefix, f'{pdb_code}_{chain_1}_{nt_number_1}_{chain_2}_{nt_number_2}_{hg_img_suffix}_90.png')
        wc_image_90 = os.path.join(base_path, wc_prefix, f'{pdb_code}_{chain_1}_{nt_number_1}_{chain_2}_{nt_number_2}_{wc_img_suffix}_90.png')

        # Extract metrics using swapped prefixes
        rscc_wc_raw = safe_float(row.get(f'RSCC_{wc_prefix}'))
        rscc_hg_raw = safe_float(row.get(f'RSCC_{hg_prefix}'))
        
        edia_wc_raw = safe_float(row.get(f'{wc_prefix}_edia'))
        edia_hg_raw = safe_float(row.get(f'{hg_prefix}_edia'))
        
        clash_bp_wc_raw = safe_float(row.get(f'{wc_prefix}_clashscore_bp'))
        clash_bp_hg_raw = safe_float(row.get(f'{hg_prefix}_clashscore_bp'))
        
        clash_neigh_wc_raw = safe_float(row.get(f'{wc_prefix}_clashscore_neighbour'))
        clash_neigh_hg_raw = safe_float(row.get(f'{hg_prefix}_clashscore_neighbour'))
        
        clash_global_wc_raw = safe_float(row.get(f'{wc_prefix}_clashscore_global'))
        clash_global_hg_raw = safe_float(row.get(f'{hg_prefix}_clashscore_global'))
        
        rwork_wc_raw = safe_float(row.get(f'r_work_{wc_prefix}'))
        rwork_hg_raw = safe_float(row.get(f'r_work_{hg_prefix}'))
        
        rfree_wc_raw = safe_float(row.get(f'r_free_{wc_prefix}'))
        rfree_hg_raw = safe_float(row.get(f'r_free_{hg_prefix}'))
        
        b_wc_raw = safe_float(row.get(f'mean_B_{wc_prefix}'))
        b_hg_raw = safe_float(row.get(f'mean_B_{hg_prefix}'))

        # Recalculate differences as WC - HG (after swap if applicable)
        rscc_wc = rscc_wc_raw
        rscc_hg = rscc_hg_raw
        rscc_diff = safe_float(rscc_wc_raw) if _is_num(rscc_wc_raw) and _is_num(rscc_hg_raw) else 'NA'
        if _is_num(rscc_diff):
            rscc_diff = float(rscc_wc_raw) - float(rscc_hg_raw)

        edia_wc = edia_wc_raw
        edia_hg = edia_hg_raw
        edia_diff = safe_float(edia_wc_raw) if _is_num(edia_wc_raw) and _is_num(edia_hg_raw) else 'NA'
        if _is_num(edia_diff):
            edia_diff = float(edia_wc_raw) - float(edia_hg_raw)

        clash_bp_wc = clash_bp_wc_raw
        clash_bp_hg = clash_bp_hg_raw
        clash_bp_diff = safe_float(clash_bp_wc_raw) if _is_num(clash_bp_wc_raw) and _is_num(clash_bp_hg_raw) else 'NA'
        if _is_num(clash_bp_diff):
            clash_bp_diff = float(clash_bp_wc_raw) - float(clash_bp_hg_raw)

        clash_neigh_wc = clash_neigh_wc_raw
        clash_neigh_hg = clash_neigh_hg_raw
        clash_neigh_diff = safe_float(clash_neigh_wc_raw) if _is_num(clash_neigh_wc_raw) and _is_num(clash_neigh_hg_raw) else 'NA'
        if _is_num(clash_neigh_diff):
            clash_neigh_diff = float(clash_neigh_wc_raw) - float(clash_neigh_hg_raw)

        clash_global_wc = clash_global_wc_raw
        clash_global_hg = clash_global_hg_raw

        rwork_wc = rwork_wc_raw
        rwork_hg = rwork_hg_raw
        
        rfree_wc = rfree_wc_raw
        rfree_hg = rfree_hg_raw

        b_wc = b_wc_raw
        b_hg = b_hg_raw

        # Calculate R-value diffs if possible
        rwork_diff = 'NA'
        rfree_diff = 'NA'
        if _is_num(rwork_wc) and _is_num(rwork_hg):
            rwork_diff = float(rwork_wc) - float(rwork_hg)
        if _is_num(rfree_wc) and _is_num(rfree_hg):
            rfree_diff = float(rfree_wc) - float(rfree_hg)

        # Check if all metrics are missing
        all_metrics = [
            rscc_wc, rscc_hg, edia_wc, edia_hg,
            clash_bp_wc, clash_bp_hg, clash_neigh_wc, clash_neigh_hg,
            clash_global_wc, clash_global_hg,
            rwork_wc, rwork_hg, rfree_wc, rfree_hg,
            b_wc, b_hg
        ]
        all_missing = all(not _is_num(x) for x in all_metrics)

        # Interpretations
        interpretations = {
            "RSCC": interpret_metric(rscc_diff, "RSCC") if _is_num(rscc_diff) else None,
            "EDIA": interpret_metric(edia_diff, "EDIA") if _is_num(edia_diff) else None,
            "Clashscore bp": interpret_metric(clash_bp_diff, "Clashscore") if _is_num(clash_bp_diff) else None,
            "Clashscore neighbour": interpret_metric(clash_neigh_diff, "Clashscore") if _is_num(clash_neigh_diff) else None,
            "R-work": interpret_metric(rwork_diff, "R-work") if _is_num(rwork_diff) else None,
            "R-free": interpret_metric(rfree_diff, "R-free") if _is_num(rfree_diff) else None,
        }

        # Core metrics for overall decision
        interpretations_core = {
            k: v for k, v in interpretations.items()
            if k in ["RSCC", "EDIA", "Clashscore bp", "Clashscore neighbour"]
        }
        overall_result = determine_overall(interpretations_core)

        # Poor density override
        poor_density = (
            _is_num(edia_hg) and _is_num(edia_wc) and (
                (float(edia_hg) < 0.5 and float(edia_wc) < 0.5)
            )
        )

        if poor_density:
            overall_result = "Poor electron density"

        # Error override - check occupancy first
        if all_missing:
            # Check occupancy file
            occupancy_file = os.path.join(PATH_TO_IMAGES, f'{pdb_code}_{chain_purine}_{nt_purine}', 'WC', 'occupancy')
            if os.path.exists(occupancy_file):
                try:
                    with open(occupancy_file, 'r') as f:
                        occupancy_value = f.read().strip()
                        # Check if occupancy is not equal to 1
                        try:
                            occ_float = float(occupancy_value)
                            if abs(occ_float - 1.0) > 0.001:
                                overall_result = "Occupancy_not_1"
                            else:
                                overall_result = "Error"
                        except ValueError:
                            overall_result = "Error"
                except Exception as e:
                    print(f"  Warning: Could not read occupancy file: {e}")
                    overall_result = "Error"
            else:
                overall_result = "Error"

        # Create a new file with metrics + classification
        RESULTS_FILE = "classification_files/classification_results.txt"
        
        # Prepare the output line with all metrics and classification
        output_line = (
            f"{pdb_code} {resolution} {chain_1} {nt_type_1} {nt_number_1} "
            f"{chain_2} {nt_type_2} {nt_number_2} "
            f"{format_value(rscc_wc)} {format_value(rscc_hg)} {format_value(rscc_diff)} "
            f"{format_value(edia_wc)} {format_value(edia_hg)} {format_value(edia_diff)} "
            f"{format_value(clash_bp_wc, 1)} {format_value(clash_bp_hg, 1)} {format_value(clash_bp_diff, 1)} "
            f"{format_value(clash_neigh_wc, 1)} {format_value(clash_neigh_hg, 1)} {format_value(clash_neigh_diff, 1)} "
            f"{format_value(clash_global_wc, 1)} {format_value(clash_global_hg, 1)} "
            f"{format_value(rwork_wc)} {format_value(rwork_hg)} {format_value(rwork_diff)} "
            f"{format_value(rfree_wc)} {format_value(rfree_hg)} {format_value(rfree_diff)} "
            f"{format_value(b_wc, 1)} {format_value(b_hg, 1)} "
            f"{conformation_bp} "
            f"{overall_result}\n"
        )
        
        # Create header if file doesn't exist
        os.makedirs("classification_files", exist_ok=True)
        if not os.path.exists(RESULTS_FILE):
            header = (
                "pdb_id resolution chain_1 nt_type_1 nt_number_1 chain_2 nt_type_2 nt_number_2 "
                "RSCC_WC RSCC_HG delta_RSCC "
                "EDIA_WC EDIA_HG delta_EDIA "
                "Clashscore_bp_WC Clashscore_bp_HG delta_Clashscore_bp "
                "Clashscore_neighbour_WC Clashscore_neighbour_HG delta_Clashscore_neighbour "
                "Clashscore_global_WC Clashscore_global_HG "
                "R_work_WC R_work_HG delta_R_work "
                "R_free_WC R_free_HG delta_R_free "
                "B_factor_WC B_factor_HG "
                "Initial_conformation "
                "classification\n"
            )
            with open(RESULTS_FILE, 'w') as f:
                f.write(header)
        
        # Append the new result
        with open(RESULTS_FILE, 'a') as f:
            f.write(output_line)

        # Create plot with 2 rows and 2 columns
        with PdfPages(OUTPUT_PDF) as pdf:
            fig = plt.figure(figsize=(14, 16))
            
            # Create grid: 2 rows for images, space at bottom for text
            gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.5], hspace=0.15, wspace=0.1)
            
            # Top row - 90 degree rotated views
            ax_hg_90 = fig.add_subplot(gs[0, 0])
            ax_wc_90 = fig.add_subplot(gs[0, 1])
            
            # Bottom row - original views
            ax_hg = fig.add_subplot(gs[1, 0])
            ax_wc = fig.add_subplot(gs[1, 1])

            # HG 90° image
            if os.path.exists(hg_image_90):
                ax_hg_90.imshow(plt.imread(hg_image_90))
                ax_hg_90.set_title(f'{purine}(syn) Model - HG (90° rotation)')
            else:
                ax_hg_90.text(0.5, 0.5, f'{purine}(syn) 90° image not found', 
                             ha='center', va='center')
            ax_hg_90.axis('off')

            # WC 90° image
            if os.path.exists(wc_image_90):
                ax_wc_90.imshow(plt.imread(wc_image_90))
                ax_wc_90.set_title(f'{purine}(anti) Model - WC (90° rotation)')
            else:
                ax_wc_90.text(0.5, 0.5, f'{purine}(anti) 90° image not found', 
                             ha='center', va='center')
            ax_wc_90.axis('off')

            # HG image
            if os.path.exists(hg_image):
                ax_hg.imshow(plt.imread(hg_image))
                ax_hg.set_title(f'{purine}(syn) Model - HG')
            else:
                ax_hg.text(0.5, 0.5, f'{purine}(syn) image not found', 
                          ha='center', va='center')
            ax_hg.axis('off')

            # WC image
            if os.path.exists(wc_image):
                ax_wc.imshow(plt.imread(wc_image))
                ax_wc.set_title(f'{purine}(anti) Model - WC')
            else:
                ax_wc.text(0.5, 0.5, f'{purine}(anti) image not found', 
                          ha='center', va='center')
            ax_wc.axis('off')

            # Format metrics text
            metric_text = f"""
PDB ID: {pdb_code} | Resolution: {resolution} Å
Chains: {chain_1}, {chain_2} | Nucleotides: {nt_type_1}.{nt_number_1} : {nt_type_2}.{nt_number_2}
Original conformation: {purine}({conformation}) - {conformation_bp}

Metric                      HG        WC        Δ(WC - HG)  Favours
------------------------------------------------------------------
RSCC                        {format_value(rscc_hg):>9} {format_value(rscc_wc):>9} {format_value(rscc_diff):>10} {interpretations['RSCC'] or 'N/A'}
EDIA                        {format_value(edia_hg):>9} {format_value(edia_wc):>9} {format_value(edia_diff):>10} {interpretations['EDIA'] or 'N/A'}
Clashscore bp               {format_value(clash_bp_hg, 1):>9} {format_value(clash_bp_wc, 1):>9} {format_value(clash_bp_diff, 1):>10} {interpretations['Clashscore bp'] or 'N/A'}
Clashscore neighbour        {format_value(clash_neigh_hg, 1):>9} {format_value(clash_neigh_wc, 1):>9} {format_value(clash_neigh_diff, 1):>10} {interpretations['Clashscore neighbour'] or 'N/A'}
Clashscore global           {format_value(clash_global_hg, 1):>9} {format_value(clash_global_wc, 1):>9}
R-work                      {format_value(rwork_hg):>9} {format_value(rwork_wc):>9} {format_value(rwork_diff):>10} {interpretations['R-work'] or 'N/A'}
R-free                      {format_value(rfree_hg):>9} {format_value(rfree_wc):>9} {format_value(rfree_diff):>10} {interpretations['R-free'] or 'N/A'}
B-factor (purine avg)       {format_value(b_hg, 1):>9} {format_value(b_wc, 1):>9}
"""
            fig.text(0.05, 0.12, metric_text, fontsize=10, family='monospace')
            fig.text(0.05, 0.08, f"Overall Result: {overall_result}", 
                    fontsize=12, family='monospace', weight='bold')

            pdf.savefig(fig)
            plt.close(fig)

        # AFTER generating the PDF, check if we need to delete the folder
        if overall_result in ["Poor electron density", "WC"]:
            folder_to_delete = base_path
            delete_folder_safely(folder_to_delete)

    except Exception as e:
        print(f" Error processing entry: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)