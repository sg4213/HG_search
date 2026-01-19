# Hoogsteen Base Pair Finder (Hoog_finder_by_bp)

A pipeline to identify and classify Watson-Crick (WC) vs Hoogsteen (HG) base pairs in X-ray crystallographic structures by comparing electron density fit metrics between both conformations.

This pipeline processes PDB structures to:
1. Download/copy PDB and MTZ files
2. Generate WC and HG conformations by flipping the purine (A or G)
3. Refine both models against the electron density
4. Calculate quality metrics (RSCC, EDIA, clashscore, R-values, B-factors)
5. Compare metrics to classify the base pair as WC, HG, or Ambiguous
6. Generate PDF reports with visualizations

## Requirements

### Software Dependencies
- **Phenix** - crystallographic refinement suite
  - `phenix.refine`
  - `phenix.clashscore`
  - `phenix.mtz2map`
- **PyMOL** - molecular visualization (with Python API)
- **EDIA scorer** - electron density assessment tool
- **Python 3** with packages:
  - `pandas`
  - `numpy`
  - `matplotlib`
- **wget** - for downloading MTZ files from PDB-REDO

### Environment Variables
export EDIA_BIN=/path/to/ediascorer      # Path to EDIA binary
export EDIA_LICENSE=<license_key>         # EDIA license key


## Input

### PairTable_X_ray.csv
A CSV file containing base pairs to analyze. Columns:

pdb_code,assembly,resolution,chi_1,chi_2,chain_1,nt_type_1,nt_number_1,chain_2,nt_type_2,nt_number_2,xxxx


### PDB Files
Pre-processed PDB files should be available at the path specified in `PDB_PATH` variable in `batch_run.sh`:
PDB_PATH="/path/to/pdb_files/"

## Usage

### Before Running (Important!)
**Clean up previous runs before starting a new analysis:**

# Remove all files in classification_files (pipeline appends to these files)
rm -f classification_files/*

# Remove folders in PDB_without_nt (existing folders with Phenix files will be skipped, not overwritten)
rm -rf PDB_without_nt/*


### Running the Full Pipeline
./batch_run.sh


This will process all entries in `PairTable_X_ray.csv`.

### Configuration
Edit `batch_run.sh` to set:

PDB_PATH="/path/to/your/pdb/files/"   # Local PDB file directory
MTZ_URL="https://pdb-redo.eu/db"       # MTZ download URL
csv_file="PairTable_X_ray.csv"         # Input CSV file


### Classification Results
`classification_files/classification_results.txt` contains the final classification with all metrics:
- RSCC (WC, HG, delta)
- EDIA (WC, HG, delta)
- Clashscore (base pair, neighbour, global)
- R-values (R-work, R-free)
- B-factors
- Initial conformation
- Final classification

| Metric | WC Favoured | HG Favoured |
|--------|-------------|-------------|
| RSCC | Δ > 0.007 | Δ < -0.007 |
| EDIA | Δ > 0.010 | Δ < -0.010 |
| Clashscore | WC < HG | HG < WC |

Overall result is based on the majority vote from RSCC, EDIA, clashscore (bp), clashscore (neighbour)

**Special cases**:
- "Poor electron density" if both EDIA values < 0.5
- "Error" if all metrics are missing
- "Occupancy_not_1" if nucleotide occupancy ≠ 1

## Individual Scripts

| Script | Description |
|--------|-------------|
| `batch_run.sh` | Main pipeline orchestrator |
| `read_PDB_MTZ_NT_AT_GT_AC.sh` | Process AT/TA/GT/TG/AC/CA/AU/UA/GU/UG base pairs |
| `read_PDB_MTZ_NT_GC.sh` | Process GC/CG base pairs |
| `flip.py` | Flip purine to HG conformation (positive residue numbers) |
| `flip_negative.py` | Flip purine to HG conformation (negative residue numbers) |
| `check_occupancy.py` | Check nucleotide occupancy |
| `protonate.py` | Add protons to structure |
| `get_rval.sh` | Calculate R-values |
| `get_rscc.sh` | Calculate RSCC |
| `get_EDIA.sh` | Calculate EDIA scores |
| `Clashes.py` | Calculate clashscores |
| `Bfactor.py` | Calculate B-factors |
| `combine_metrics.py` | Merge all metrics |
| `make_report.py` | Generate classification and PDF report |

