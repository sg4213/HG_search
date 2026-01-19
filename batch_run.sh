#!/bin/bash
set -euo pipefail

#PDB_URL="https://pdb-redo.eu/db"
#PDB_PATH="/mnt/hdd_04/ec3867/NAFinder/NAFinder_05122025/X-ray/pdb"
PDB_PATH="/mnt/hdd_04/ec3867/NAFinder/NAFinder_20260108/X-ray/pdb_dssr/"
MTZ_URL="https://pdb-redo.eu/db"
csv_file="PairTable_X_ray.csv"   # or pass as $1 if you prefer

#download_pdb() {
#  local pdbid="$1"
#  local url="${PDB_URL}/${pdbid}/${pdbid}_final.pdb"
#  local out="${pdbid}_final.pdb"
#  echo "Downloading PDB: ${url} -> ${out}"
#  wget -q -c -O "${out}" "${url}" || echo "Failed to download ${url}"
#}

download_mtz() {
  local pdbid="$1"
  local url="${MTZ_URL}/${pdbid}/${pdbid}_final.mtz"
  local out="${pdbid}_final.mtz"
  #echo "Downloading MTZ: ${url} -> ${out}"
  wget -q -c -O "${out}" "${url}" || echo "Failed to download ${url}"
}

# CSV columns:
# pdb_code, assembly, chain_1, nt_type_1, nt_number_1, chain_2, nt_type_2, nt_number_2, xxxx
while IFS=, read -r pdb_code assembly reso chi_1 chi_2 chain_1 nt_type_1 nt_number_1 chain_2 nt_type_2 nt_number_2 xxxx; do
  [[ -z "${pdb_code:-}" ]] && continue
  if [[ "${pdb_code}" == "pdb_code" ]]; then
    continue
  fi

  #download_pdb "${pdb_code}"
  cp "$PDB_PATH/${pdb_code}.pdb${assembly}" "."
  mv "${pdb_code}.pdb${assembly}" "${pdb_code}_final.pdb" 
  download_mtz "${pdb_code}"

  pair="${nt_type_1}${nt_type_2}"
  script=""
  case "${pair}" in
    AT|TA|TG|GT|AC|CA|AU|UA|GU|UG|GC|CG)
      script="./read_PDB_MTZ_NT_AT_GT_AC.sh"
      ;;
    GC|CG)
      script="./read_PDB_MTZ_NT_GC.sh"
      ;;
    *)
      echo "Other bp: ${pair} — no script mapped"
      ;;
  esac

  if [[ -n "${script}" ]]; then
    cmd="${script} ${pdb_code} ${pdb_code} ${chain_1} ${nt_type_1} ${nt_number_1} ${chain_2} ${nt_type_2} ${nt_number_2} ${xxxx}"
    echo "Running: ${cmd}"
    if ${cmd}; then
      echo "OK: ${cmd}" >> out.txt
    else
      echo "ERROR: ${cmd}" >> out_error.txt
    fi
  fi
  echo "Rval"
  ./get_rval.sh ${pdb_code} ${chain_1} ${nt_type_1} ${nt_number_1} ${chain_2} ${nt_type_2} ${nt_number_2}
  echo "RSCC"
  ./get_rscc.sh ${pdb_code} ${chain_1} ${nt_type_1} ${nt_number_1} ${chain_2} ${nt_type_2} ${nt_number_2}
  echo "EDIA"
  ./get_EDIA.sh ${pdb_code} ${chain_1} ${nt_type_1} ${nt_number_1} ${chain_2} ${nt_type_2} ${nt_number_2}
  echo "Clashscore"
  python3 Clashes.py ${pdb_code} ${chain_1} ${nt_type_1} ${nt_number_1} ${chain_2} ${nt_type_2} ${nt_number_2}
  echo "B-facor"
  python3 Bfactor.py ${pdb_code} ${chain_1} ${nt_type_1} ${nt_number_1} ${chain_2} ${nt_type_2} ${nt_number_2}
  echo "Combine metrics"
  python3 combine_metrics.py ${pdb_code} ${chain_1} ${nt_type_1} ${nt_number_1} ${chain_2} ${nt_type_2} ${nt_number_2}
  echo "Make report"
  python3 make_report.py ${pdb_code} ${reso} ${chain_1} ${nt_type_1} ${nt_number_1} ${chain_2} ${nt_type_2} ${nt_number_2} ${chi_1} ${chi_2} 
  rm "${pdb_code}_final.pdb"
  rm "${pdb_code}_final.mtz"

done < "${csv_file}"
