#! /bin/bash
mkdir Bad_Data
file="b_factor.txt"
  while IFS='	' read -r f1 f2
  do
    echo $f2
    if [[ "$f2" -gt "70" ]]
    then
      echo "Bad data"
      mv PDB_without_nt/${f1::4}_${f1:7:9} Bad_Data
    fi
    echo "Good data"

  done < $file
