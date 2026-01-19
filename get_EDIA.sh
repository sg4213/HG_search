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

EDIA_BIN="${EDIA_BIN:-/home/sg4109/software/EDIA/ediascorer_1.1.0/ediascorer}"
EDIA_LICENSE="${EDIA_LICENSE:-AAAAAAAAljfUAAAAU2CEpBAtTlNi83vbDe9jtzbdHCo8=}"
BASE="PDB_without_nt"

# Dependency check
if ! command -v "$EDIA_BIN" >/dev/null 2>&1; then
  echo "ERROR: ediascorer not found at: $EDIA_BIN (override with EDIA_BIN=...)" >&2
  exit 1
fi

pdbid="$pdb_id"
chain="$chain_purine"
nt="$nt_purine"

# Locate tuple directory
target="$(find "$BASE" -type d -name "${pdbid}_${chain}_${nt}" | head -n1 || true)"
if [[ -z "$target" ]]; then
  echo "ERROR: directory ${pdbid}_${chain}_${nt} not found under $BASE" >&2
  exit 1
fi

OUTFILE="${target}/ediascorer_log.txt"
: > "$OUTFILE"  # truncate

log() { printf '%s %s\n' "$(date '+%F %T')" "$*" | tee -a "$OUTFILE" >/dev/null; }

run_edia() {
  local dir="$1"
  local pdb_glob="$2"
  local ccp4_glob="$3"
  
  if [[ ! -d "$dir" ]]; then
    log "WARNING: ${dir} not found — skipping"
    return 0
  fi
  
  local outdir="${dir}/edia_out"
  
  # Check if output already exists
  if [[ -d "$outdir" ]]; then
    local score_file
    score_file="$(compgen -G "${outdir}/*_001structurescores.csv" | head -n1 || true)"
    if [[ -n "$score_file" ]]; then
      log "SKIPPING: EDIA output already exists in ${outdir}"
      return 0
    fi
  fi
  
  # Pick first match for each
  local pdb ccp4
  pdb="$(compgen -G "${dir}/${pdb_glob}" | head -n1 || true)"
  ccp4="$(compgen -G "${dir}/${ccp4_glob}" | head -n1 || true)"
  
  if [[ -z "$pdb" || -z "$ccp4" ]]; then
    log "ERROR: missing PDB or CCP4 in ${dir}"
    return 0
  fi
  
  mkdir -p "$outdir"
  chmod 777 -R "$outdir" 2>/dev/null || true
  
  log "Running EDIA in ${dir}"
  if "$EDIA_BIN" \
        -l "$EDIA_LICENSE" \
        -t "$pdb" \
        -d "$ccp4" \
        -o "${outdir}/" >>"$OUTFILE" 2>&1; then
    log "SUCCESS: EDIA completed for ${dir}"
  else
    log "ERROR: EDIA failed for ${dir}"
  fi
}

log "Target: ${target}"

# ---- WC ----
run_edia "${target}/WC" "*final_refine_001.pdb" "*final_refine_001_2mFo-DFc.ccp4"

# ---- HG ----
run_edia "${target}/HG" "*final_flipped_refine_001_refine_001.pdb" "*final_flipped_refine_001_refine_001_2mFo-DFc.ccp4"

# ---------------------------
# extract_nt_edia:
# find line "r,*,<nt_purine>,<chain_purine>,..." and take column 5
# ---------------------------
extract_nt_edia() {
  local parent_dir="$1"
  local __resultvar="$2"
  local outdir="${parent_dir}/edia_out"
  
  if [[ ! -d "$outdir" ]]; then
    log "WARNING: ${outdir} not found — skipping EDIA extraction"
    printf -v "$__resultvar" 'None'
    return 0
  fi
  
  local score_file
  score_file="$(compgen -G "${outdir}/*_001structurescores.csv" | head -n1 || true)"
  
  if [[ -z "$score_file" ]]; then
    log "WARNING: no EDIA score file found in ${outdir}"
    printf -v "$__resultvar" 'None'
    return 0
  fi
  
  local value
  value="$(
    awk -F',' -v nt="$nt_purine" -v ch="$chain_purine" '
      ($3 == nt && $4 == ch) { print $5; exit }
    ' "$score_file"
  )"
  
  if [[ -z "$value" ]]; then
    log "WARNING: could not find EDIA entry for nt=${nt_purine}, chain=${chain_purine} in ${score_file}"
    printf -v "$__resultvar" 'None'
    return 0
  fi
  
  printf -v "$__resultvar" '%s' "$value"
}

# ---------------------------
# Get WC_edia and HG_edia
# ---------------------------
WC_edia="None"
HG_edia="None"

extract_nt_edia "${target}/WC" WC_edia
extract_nt_edia "${target}/HG" HG_edia

summary_file="classification_files/EDIA_summary.txt"

if [[ ! -f "$summary_file" ]]; then
  echo "pdb_id chain_1 nt_type_1 nt_number_1 chain_2 nt_type_2 nt_number_2 WC_edia HG_edia" \
    > "$summary_file"
fi

echo "$pdb_id $chain_1 $nt_type_1 $nt_number_1 $chain_2 $nt_type_2 $nt_number_2 $WC_edia $HG_edia" \
  >> "$summary_file"