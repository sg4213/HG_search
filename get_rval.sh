#!/usr/bin/env bash
set -euo pipefail

if [ $# -ne 7 ]; then
    echo "Usage: $0 pdb_id chain_1 nt_type_1 nt_number_1 chain_2 nt_type_2 nt_number_2" >&2
    exit 1
fi

pdb_id="$1"
#pdb_file="$1_final.pdb"
mtz_file="$1_final.mtz"
chain_1="$2"
nt_type_1="$3"
nt_number_1="$4"
chain_2="$5"
nt_type_2="$6"
nt_number_2="$7"

# ---------------------------
# Choose purine (A/G)
# ---------------------------
if [[ "$nt_type_1" == "G" || "$nt_type_1" == "A" ]]; then
    chain_purine="$chain_1"
    nt_purine="$nt_number_1"
elif [[ "$nt_type_2" == "G" || "$nt_type_2" == "A" ]]; then
    chain_purine="$chain_2"
    nt_purine="$nt_number_2"
else
    # mismatch or modified, nothing to do
    exit 1
fi

# ---------------------------
# Find target directory
# ---------------------------
base="PDB_without_nt"
target=$(find "$base" -type d -name "${pdb_id}_${chain_purine}_${nt_purine}" | head -n1)

if [[ -z "$target" ]]; then
    exit 1
fi

cd "$target"

# ---------------------------
# Extract R-values
# ---------------------------
extract_r_values() {
  local file=$1
  awk '
    /R VALUE.*WORKING.*TEST/        { r_total = $NF }
    /R VALUE.*WORKING SET/          { r_work  = $NF }
    # match only the lone “FREE R VALUE” (no other words after)
    /^REMARK.*FREE R VALUE[[:space:]]*:[[:space:]]*[0-9.]+$/ {
                                     r_free  = $NF }
    END {
      printf "%s %s %s\n", r_total, r_work, r_free
    }
  ' "$file"
}


summary_file="../../classification_files/R_values_summary.txt"

if [[ ! -f "$summary_file" ]]; then
  echo "pdb_id chain_1 nt_type_1 nt_number_1 chain_2 nt_type_2 nt_number_2 r_total_WC r_work_WC r_free_WC r_total_HG r_work_HG r_free_HG" \
    > "$summary_file"
fi

# ---------------------------
# Collect WC / HG R-values
# ---------------------------
r_total_WC="NA"; r_work_WC="NA"; r_free_WC="NA"
r_total_HG="NA"; r_work_HG="NA"; r_free_HG="NA"

for mode in WC HG; do

  cd "$mode"

  if [[ "$mode" == "WC" ]]; then
    pdbfile="${pdb_id}_final_refine_001.pdb"
    if [[ -f "$pdbfile" ]]; then
      read r_total_WC r_work_WC r_free_WC < <(extract_r_values "$pdbfile")
    else
      echo "WARNING: WC file '$pdbfile' not found" >&2
    fi
  else  # HG
    pdbfile="${pdb_id}_final_flipped_refine_001_refine_001.pdb"
    if [[ -f "$pdbfile" ]]; then
      read r_total_HG r_work_HG r_free_HG < <(extract_r_values "$pdbfile")
    else
      echo "WARNING: HG file '$pdbfile' not found" >&2
    fi
  fi

  cd ..
done

# ---------------------------
# Append one line to summary
# ---------------------------
echo "$pdb_id $chain_1 $nt_type_1 $nt_number_1 $chain_2 $nt_type_2 $nt_number_2 $r_total_WC $r_work_WC $r_free_WC $r_total_HG $r_work_HG $r_free_HG" \
  >> "$summary_file"
