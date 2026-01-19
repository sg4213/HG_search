#!/usr/bin/env python3
import os
import sys
import pandas as pd
import numpy as np


SUMMARY_FILES = {
    'bfactor': 'classification_files/Bfactor_summary.txt',
    'clashscore': 'classification_files/clashscore_summary.txt',
    'edia': 'classification_files/EDIA_summary.txt',
    'rscc': 'classification_files/RSCC_summary.txt',
    'rvalues': 'classification_files/R_values_summary.txt',
}

OUT_COMBINED = "classification_files/combined_metrics.txt"

# Expected column order
EXPECTED_COLUMNS = [
    'pdb_id', 'chain_1', 'nt_type_1', 'nt_number_1', 'chain_2', 'nt_type_2', 'nt_number_2',
    'r_total_WC', 'r_work_WC', 'r_free_WC', 'r_total_HG', 'r_work_HG', 'r_free_HG',
    'RSCC_WC', 'RSCC_HG', 'delta_RSCC',
    'WC_edia', 'HG_edia', 'delta_edia',
    'WC_clashscore_global', 'HG_clashscore_global',
    'WC_clashscore_bp', 'HG_clashscore_bp', 'delta_clashscore_bp',
    'WC_clashscore_neighbour', 'HG_clashscore_neighbour', 'delta_clashscore_neighbour',
    'mean_B_HG', 'mean_B_WC'
]

# -------------------------
# Helper Functions
# -------------------------
def read_summary_file(filepath):
    """Read a summary file and return a DataFrame."""
    if not os.path.exists(filepath):
        return None
    
    try:
        df = pd.read_csv(filepath, sep=r'\s+')
        return df
    except Exception as e:
        print(f" Error reading {filepath}: {e}", file=sys.stderr)
        return None

def safe_subtract(val1, val2):
    """Safely subtract two values, handling None/NaN."""
    try:
        v1 = pd.to_numeric(val1, errors='coerce')
        v2 = pd.to_numeric(val2, errors='coerce')
        if pd.isna(v1) or pd.isna(v2):
            return None
        return v1 - v2
    except:
        return None

def get_value_from_df(df, pdb_id, chain_1, nt_type_1, nt_number_1, chain_2, nt_type_2, nt_number_2, column):
    """Get a specific value from a dataframe based on the key columns."""
    if df is None or column not in df.columns:
        return None
    
    mask = (
        (df['pdb_id'] == pdb_id) &
        (df['chain_1'] == chain_1) &
        (df['nt_type_1'] == nt_type_1) &
        (df['nt_number_1'] == nt_number_1) &
        (df['chain_2'] == chain_2) &
        (df['nt_type_2'] == nt_type_2) &
        (df['nt_number_2'] == nt_number_2)
    )
    
    matches = df.loc[mask, column]
    if len(matches) > 0:
        val = matches.iloc[0]
        return None if pd.isna(val) else val
    return None

def create_combined_table():
    """Create the combined table from all summary files."""
    
    # Read all summary files
    df_bfactor = read_summary_file(SUMMARY_FILES['bfactor'])
    df_clashscore = read_summary_file(SUMMARY_FILES['clashscore'])
    df_edia = read_summary_file(SUMMARY_FILES['edia'])
    df_rscc = read_summary_file(SUMMARY_FILES['rscc'])
    df_rvalues = read_summary_file(SUMMARY_FILES['rvalues'])
    
    # Check if at least one file exists
    if all(df is None for df in [df_bfactor, df_clashscore, df_edia, df_rscc, df_rvalues]):
        print(" Error: No summary files found", file=sys.stderr)
        return False
    
    # Start with the first available dataframe as base
    base_df = None
    for df in [df_rvalues, df_rscc, df_edia, df_clashscore, df_bfactor]:
        if df is not None:
            base_df = df[['pdb_id', 'chain_1', 'nt_type_1', 'nt_number_1', 
                          'chain_2', 'nt_type_2', 'nt_number_2']].copy()
            break
    
    if base_df is None:
        print(" Error: Could not establish base dataframe", file=sys.stderr)
        return False
    
    # Initialize result dataframe with base columns
    result_df = base_df.copy()
    
    # Add R-values columns
    if df_rvalues is not None:
        for col in ['r_total_WC', 'r_work_WC', 'r_free_WC', 'r_total_HG', 'r_work_HG', 'r_free_HG']:
            if col in df_rvalues.columns:
                result_df = result_df.merge(
                    df_rvalues[['pdb_id', 'chain_1', 'nt_type_1', 'nt_number_1', 
                               'chain_2', 'nt_type_2', 'nt_number_2', col]],
                    on=['pdb_id', 'chain_1', 'nt_type_1', 'nt_number_1', 
                        'chain_2', 'nt_type_2', 'nt_number_2'],
                    how='left'
                )
    
    # Add RSCC columns
    if df_rscc is not None:
        for col in ['RSCC_WC', 'RSCC_HG']:
            if col in df_rscc.columns:
                result_df = result_df.merge(
                    df_rscc[['pdb_id', 'chain_1', 'nt_type_1', 'nt_number_1', 
                            'chain_2', 'nt_type_2', 'nt_number_2', col]],
                    on=['pdb_id', 'chain_1', 'nt_type_1', 'nt_number_1', 
                        'chain_2', 'nt_type_2', 'nt_number_2'],
                    how='left'
                )
    
    # Add EDIA columns
    if df_edia is not None:
        for col in ['WC_edia', 'HG_edia']:
            if col in df_edia.columns:
                result_df = result_df.merge(
                    df_edia[['pdb_id', 'chain_1', 'nt_type_1', 'nt_number_1', 
                            'chain_2', 'nt_type_2', 'nt_number_2', col]],
                    on=['pdb_id', 'chain_1', 'nt_type_1', 'nt_number_1', 
                        'chain_2', 'nt_type_2', 'nt_number_2'],
                    how='left'
                )
    
    # Add Clashscore columns
    if df_clashscore is not None:
        for col in ['WC_clashscore_global', 'HG_clashscore_global', 
                    'WC_clashscore_bp', 'HG_clashscore_bp',
                    'WC_clashscore_neighbour', 'HG_clashscore_neighbour']:
            if col in df_clashscore.columns:
                result_df = result_df.merge(
                    df_clashscore[['pdb_id', 'chain_1', 'nt_type_1', 'nt_number_1', 
                                  'chain_2', 'nt_type_2', 'nt_number_2', col]],
                    on=['pdb_id', 'chain_1', 'nt_type_1', 'nt_number_1', 
                        'chain_2', 'nt_type_2', 'nt_number_2'],
                    how='left'
                )
    
    # Add B-factor columns
    if df_bfactor is not None:
        for col in ['mean_B_HG', 'mean_B_WC']:
            if col in df_bfactor.columns:
                result_df = result_df.merge(
                    df_bfactor[['pdb_id', 'chain_1', 'nt_type_1', 'nt_number_1', 
                               'chain_2', 'nt_type_2', 'nt_number_2', col]],
                    on=['pdb_id', 'chain_1', 'nt_type_1', 'nt_number_1', 
                        'chain_2', 'nt_type_2', 'nt_number_2'],
                    how='left'
                )
    
    # Remove duplicate rows
    result_df = result_df.drop_duplicates(
        subset=['pdb_id', 'chain_1', 'nt_type_1', 'nt_number_1', 
                'chain_2', 'nt_type_2', 'nt_number_2']
    )
    
    # Calculate difference columns
    if 'RSCC_WC' in result_df.columns and 'RSCC_HG' in result_df.columns:
        result_df['delta_RSCC'] = result_df.apply(
            lambda row: safe_subtract(row['RSCC_WC'], row['RSCC_HG']), axis=1
        )
    
    if 'WC_edia' in result_df.columns and 'HG_edia' in result_df.columns:
        result_df['delta_edia'] = result_df.apply(
            lambda row: safe_subtract(row['WC_edia'], row['HG_edia']), axis=1
        )
    
    if 'WC_clashscore_bp' in result_df.columns and 'HG_clashscore_bp' in result_df.columns:
        result_df['delta_clashscore_bp'] = result_df.apply(
            lambda row: safe_subtract(row['WC_clashscore_bp'], row['HG_clashscore_bp']), axis=1
        )
    
    if 'WC_clashscore_neighbour' in result_df.columns and 'HG_clashscore_neighbour' in result_df.columns:
        result_df['delta_clashscore_neighbour'] = result_df.apply(
            lambda row: safe_subtract(row['WC_clashscore_neighbour'], row['HG_clashscore_neighbour']), axis=1
        )
    
    # Ensure all expected columns exist
    for col in EXPECTED_COLUMNS:
        if col not in result_df.columns:
            result_df[col] = None
    
    # Reorder columns
    result_df = result_df[EXPECTED_COLUMNS]
    
    # Write to output file
    result_df.to_csv(OUT_COMBINED, sep=' ', index=False, na_rep='None')
    
    return True

def add_or_update_entry(pdb_id, chain_1, nt_type_1, nt_number_1, chain_2, nt_type_2, nt_number_2):
    """Add or update a single entry in the combined metrics table."""
    
    # Create output directory if needed
    os.makedirs("classification_files", exist_ok=True)
    
    # Check if combined table exists
    if not os.path.exists(OUT_COMBINED):
        print(f" Combined table not found at {OUT_COMBINED}")
        if not create_combined_table():
            print(" Failed to create combined table", file=sys.stderr)
            sys.exit(1)
    
    # Read existing combined table
    try:
        combined_df = pd.read_csv(OUT_COMBINED, sep=r'\s+')
    except Exception as e:
        print(f" Error reading combined table: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Read all summary files
    df_bfactor = read_summary_file(SUMMARY_FILES['bfactor'])
    df_clashscore = read_summary_file(SUMMARY_FILES['clashscore'])
    df_edia = read_summary_file(SUMMARY_FILES['edia'])
    df_rscc = read_summary_file(SUMMARY_FILES['rscc'])
    df_rvalues = read_summary_file(SUMMARY_FILES['rvalues'])
    
    # Build the new row
    new_row = {
        'pdb_id': pdb_id,
        'chain_1': chain_1,
        'nt_type_1': nt_type_1,
        'nt_number_1': nt_number_1,
        'chain_2': chain_2,
        'nt_type_2': nt_type_2,
        'nt_number_2': nt_number_2,
    }
    
    # Get R-values
    for col in ['r_total_WC', 'r_work_WC', 'r_free_WC', 'r_total_HG', 'r_work_HG', 'r_free_HG']:
        new_row[col] = get_value_from_df(df_rvalues, pdb_id, chain_1, nt_type_1, nt_number_1, 
                                         chain_2, nt_type_2, nt_number_2, col)
    
    # Get RSCC
    for col in ['RSCC_WC', 'RSCC_HG']:
        new_row[col] = get_value_from_df(df_rscc, pdb_id, chain_1, nt_type_1, nt_number_1, 
                                         chain_2, nt_type_2, nt_number_2, col)
    
    # Get EDIA
    for col in ['WC_edia', 'HG_edia']:
        new_row[col] = get_value_from_df(df_edia, pdb_id, chain_1, nt_type_1, nt_number_1, 
                                         chain_2, nt_type_2, nt_number_2, col)
    
    # Get Clashscores
    for col in ['WC_clashscore_global', 'HG_clashscore_global', 
                'WC_clashscore_bp', 'HG_clashscore_bp',
                'WC_clashscore_neighbour', 'HG_clashscore_neighbour']:
        new_row[col] = get_value_from_df(df_clashscore, pdb_id, chain_1, nt_type_1, nt_number_1, 
                                         chain_2, nt_type_2, nt_number_2, col)
    
    # Get B-factors
    for col in ['mean_B_HG', 'mean_B_WC']:
        new_row[col] = get_value_from_df(df_bfactor, pdb_id, chain_1, nt_type_1, nt_number_1, 
                                         chain_2, nt_type_2, nt_number_2, col)
    
    # Calculate deltas
    new_row['delta_RSCC'] = safe_subtract(new_row.get('RSCC_WC'), new_row.get('RSCC_HG'))
    new_row['delta_edia'] = safe_subtract(new_row.get('WC_edia'), new_row.get('HG_edia'))
    new_row['delta_clashscore_bp'] = safe_subtract(new_row.get('WC_clashscore_bp'), new_row.get('HG_clashscore_bp'))
    new_row['delta_clashscore_neighbour'] = safe_subtract(new_row.get('WC_clashscore_neighbour'), new_row.get('HG_clashscore_neighbour'))
    
    # Check if entry already exists
    mask = (
        (combined_df['pdb_id'] == pdb_id) &
        (combined_df['chain_1'] == chain_1) &
        (combined_df['nt_type_1'] == nt_type_1) &
        (combined_df['nt_number_1'] == nt_number_1) &
        (combined_df['chain_2'] == chain_2) &
        (combined_df['nt_type_2'] == nt_type_2) &
        (combined_df['nt_number_2'] == nt_number_2)
    )
    
    if mask.any():
        # Update existing entry
        for col, val in new_row.items():
            combined_df.loc[mask, col] = val

    else:
        # Add new entry
        new_df = pd.DataFrame([new_row])
        combined_df = pd.concat([combined_df, new_df], ignore_index=True)
    
    # Ensure all expected columns exist
    for col in EXPECTED_COLUMNS:
        if col not in combined_df.columns:
            combined_df[col] = None
    
    # Reorder columns
    combined_df = combined_df[EXPECTED_COLUMNS]
    
    # Write back to file
    combined_df.to_csv(OUT_COMBINED, sep=' ', index=False, na_rep='None')

# -------------------------
# Main
# -------------------------
if __name__ == "__main__":
    if len(sys.argv) == 1:
        # No arguments: create/rebuild entire combined table
        os.makedirs("classification_files", exist_ok=True)
        if create_combined_table():
            None
        else:
            print("\n Failed to create combined metrics table", file=sys.stderr)
            sys.exit(1)
    
    elif len(sys.argv) == 8:
        # 7 arguments: add/update single entry
        pdb_id = sys.argv[1]
        chain_1 = sys.argv[2]
        nt_type_1 = sys.argv[3]
        nt_number_1 = int(sys.argv[4])
        chain_2 = sys.argv[5]
        nt_type_2 = sys.argv[6]
        nt_number_2 = int(sys.argv[7])
        
        add_or_update_entry(pdb_id, chain_1, nt_type_1, nt_number_1, chain_2, nt_type_2, nt_number_2)
    
    else:
        sys.exit(1)