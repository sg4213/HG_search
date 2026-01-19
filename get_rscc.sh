#!/usr/bin/env bash
set -euo pipefail

if [ $# -ne 7 ]; then
    echo "Usage: $0 pdb_id chain_1 nt_type_1 nt_number_1 chain_2 nt_type_2 nt_number_2" >&2
    exit 1
fi

pdb_id="$1"
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
# Function to run RSCC in subdir
# ---------------------------
run_rscc() {
  local subdir=$1
  local pdb_glob=$2
  local mtz_glob=$3

  if [[ ! -d "$subdir" ]]; then
    return 0
  fi

  pushd "$subdir" > /dev/null

  # Pick the first matching file for each pattern
  local pdb mtz
  pdb="$(compgen -G "$pdb_glob" | head -n1 || true)"
  mtz="$(compgen -G "$mtz_glob" | head -n1 || true)"

  if [[ -z "$pdb" || -z "$mtz" ]]; then
    popd > /dev/null
    return 0
  fi

  phenix.real_space_correlation "$pdb" "$mtz" > RSCC_report.txt
  popd > /dev/null
}

# WC: *_final_refine_001.{pdb,mtz}
run_rscc "WC" "*_final_refine_001.pdb" "*_final_refine_001.mtz"

# HG: *_final_flipped_refine_001_refine_001.{pdb,mtz}
run_rscc "HG" "*_final_flipped_refine_001_refine_001.pdb" "*_final_flipped_refine_001_refine_001.mtz"


extract_nt_cc() {
  local file="$1"
  local chain="$2"
  local nt="$3"

  if [[ ! -f "$file" ]]; then
    echo "NA"
    return 0
  fi

  awk -v chain="$chain" -v nt="$nt" '
    BEGIN {
      header_type = ""
      header_seen = 0
      data_start  = 0
    }

    # detect header markers
    index($0, "<----id string---->") > 0 {
      header_type = "long"
      header_seen = 1
      next
    }

    index($0, "<id string>") > 0 && index($0, "----") == 0 {
      header_type = "short"
      header_seen = 1
      next
    }

    # the line immediately after the header marker is the column header → skip once
    header_seen && !data_start {
      data_start = 1
      next
    }

    # still before data
    !data_start { next }

    # skip empty lines
    NF == 0 { next }

    {
      chain_f = 0
      id_f    = 0
      cc_f    = 0

      if (header_type == "short") {
        # short header: either 9 or 8 columns
        #   9: chain altloc resname id occ ADP CC Rho1 Rho2
        #   8: chain resname id occ ADP CC Rho1 Rho2
        if (NF == 9) {
          chain_f = 1
          id_f    = 4
          cc_f    = 7
        } else if (NF == 8) {
          chain_f = 1
          id_f    = 3
          cc_f    = 6
        } else {
          next
        }
      } else if (header_type == "long") {
        # long header: either 10 or 9 columns
        #   10: chain altloc resname id atom occ ADP CC Rho1 Rho2
        #   9:  chain resname id atom occ ADP CC Rho1 Rho2
        if (NF == 10) {
          chain_f = 1
          id_f    = 4
          cc_f    = 8
        } else if (NF == 9) {
          chain_f = 1
          id_f    = 3
          cc_f    = 7
        } else {
          next
        }
      } else {
        # header not recognized yet
        next
      }

      # match this nucleotide
      if ($chain_f == chain && $(id_f) == nt) {
        sum += $(cc_f)
        n++
      }
    }

    END {
      if (n > 0) {
        printf "%.4f\n", sum / n
      }
    }
  ' "$file"
}


RSCC_WC="NA"
RSCC_HG="NA"

if [[ -f "WC/RSCC_report.txt" ]]; then
  RSCC_WC=$(extract_nt_cc "WC/RSCC_report.txt" "$chain_purine" "$nt_purine")
fi

if [[ -f "HG/RSCC_report.txt" ]]; then
  RSCC_HG=$(extract_nt_cc "HG/RSCC_report.txt" "$chain_purine" "$nt_purine")
fi




summary_file="../../classification_files/RSCC_summary.txt"

if [[ ! -f "$summary_file" ]]; then
  echo "pdb_id chain_1 nt_type_1 nt_number_1 chain_2 nt_type_2 nt_number_2 RSCC_WC RSCC_HG" > "$summary_file"
fi

echo "$pdb_id $chain_1 $nt_type_1 $nt_number_1 $chain_2 $nt_type_2 $nt_number_2 $RSCC_WC $RSCC_HG" >> "$summary_file"
