import sys, re
from collections import defaultdict

def split_resi(s: str):
    # Allow optional leading '-' for negative residue numbers
    m = re.fullmatch(r"(-?\d+)([A-Za-z]?)", s)
    if not m:
        sys.exit("RESI must look like -1, 0, 123 or 123A")
    return int(m.group(1)), (m.group(2) or " ")

pdb, chain, resi = sys.argv[1], sys.argv[2], sys.argv[3]
resn, icode = split_resi(resi)

rows, by_alt = [], defaultdict(list)
with open(pdb) as fh:
    for ln in fh:
        if not (ln.startswith("ATOM  ") or ln.startswith("HETATM")):
            continue
        if ln[21] != chain:                          # chain ID
            continue
        try:
            seq = int(ln[22:26])                     # resSeq (can be negative)
        except ValueError:
            continue
        ic = ln[26]                                  # insertion code

        if seq != resn or (icode != " " and ic != icode):
            continue

        resname = ln[17:20].strip()
        atom    = ln[12:16].strip()
        alt     = (ln[16].strip() or ".")           # altLoc
        occ     = float(ln[54:60])                  # occupancy
        rows.append((resname, f"{seq}{'' if ic==' ' else ic}", atom, alt, occ))
        by_alt[alt].append(occ)

if not rows:
    sys.exit("No atoms found for that selection.")

for alt, qs in by_alt.items():
    print(f"{sum(qs)/len(qs):.3f}")
