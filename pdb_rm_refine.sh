#! /bin/bash
rm -R PDB
rm -R PDB_without_nt
mkdir PDB


# Table with pdb: test_pdb.csv


PROGNAME=$0
PDB_URL="https://files.rcsb.org/download"
MTZ_URL="https://edmaps.rcsb.org/coefficients"

usage() {
  cat << EOF >&2
Usage: $PROGNAME -f <file> [-o <dir>] [-c] [-p]

 -f <file>: the input file containing a comma-separated list of PDB ids
 -o  <dir>: the output dir, default: current dir
 -p       : download a pdb.gz file for each PDB id (not available for large structures)
 -m       : download a cif.gz file for each PDB id

EOF
  exit 1
}

download_pdb() {
  url_pdb="$PDB_URL/$1"
  out=$2/$1
  echo "Downloading $url_pdb to $out"
  curl -s -f $url_pdb -o $out || echo "Failed to download $url_pdb"
}

download_mtz() {
  url_mtz="$MTZ_URL/$1"
  out=$2
  echo "Downloading $url_mtz to $out"
  wget $url_mtz -P $out|| echo "Failed to download $url_mtz"
}
select_nt() {
  if [ "$1" == "DC" ]
  then
       echo "$4" > delete.txt
       nt=${4:2:9}
  elif [ "$1" == "DG" ]
  then
       echo "$3" > delete.txt
       nt=${3:2:9}

  elif [ "$1" == "DA" ]
  then
       echo "$3" > delete.txt
       nt=${3:2:9}
  elif [ "$1" == "DT" ]
  then
       echo "$4" > delete.txt
       nt=${4:2:9}
  else    echo "mismatch or modification"
  fi
}

########################################
# Download a list of PDBs
########################################


awk -F';' '{print $1}' TATA.csv >> d.txt
tr '\n' ',' < d.txt > download.txt
rm d.txt

listfile="download.txt"
outdir="PDB"

while getopts f:o:cpaxsmr o
do
  case $o in
    (f) listfile=$OPTARG;;
    (o) outdir=$OPTARG;;
    (*) usage
  esac
done
shift "$((OPTIND - 1))"

contents=$(cat $listfile)

# see https://stackoverflow.com/questions/918886/how-do-i-split-a-string-on-a-delimiter-in-bash#tab-top
IFS=',' read -ra tokens <<< "$contents"

for token in "${tokens[@]}"
do
  download_pdb ${token}.pdb $outdir
  download_mtz ${token}.mtz $outdir
done


########################################
  # Delete nt
########################################



mkdir PDB_without_nt

file="TATA.csv"
  while IFS=';' read -r f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18
  do

    select_nt $f15 $f16 $f11 $f12
    mkdir PDB_without_nt/$f1\_$nt
    cp PDB/$f1.mtz PDB_without_nt/$f1\_$nt
    #echo PDB/$f1.pdb
    echo ${nt:4:6}
    n=${nt:4:6}
    #echo $n
    number=${#n}
    echo $number
    if [ "$number" == 1 ]
    then
      sed "/ ${nt::1}   ${nt:4:6}  /d" PDB/$f1.pdb > PDB_without_nt/$f1\_$nt/$f1\_$nt.pdb
    elif [ "$number" == 2 ]
    then
      sed "/ ${nt::1}  ${nt:4:6} /d" PDB/$f1.pdb > PDB_without_nt/$f1\_$nt/$f1\_$nt.pdb
    elif [ "$number" == 3 ]
    then
      sed "/ ${nt::1} ${nt:4:6} /d" PDB/$f1.pdb > PDB_without_nt/$f1\_$nt/$f1\_$nt.pdb
    elif [ "$number" == "4" ]
    then
      sed "/ ${nt::1}${nt:4:6} /d" PDB/$f1.pdb > PDB_without_nt/$f1\_$nt/$f1\_$nt.pdb
    fi

  done < $file


########################################
  # Run Phenix refine
########################################


 cd PDB_without_nt
 echo * >> folders.txt

 contents=$(cat folders.txt)

 # see https://stackoverflow.com/questions/918886/how-do-i-split-a-string-on-a-delimiter-in-bash#tab-top
 IFS=' ' read -ra tokens <<< "$contents"

 for token in "${tokens[@]}"
 do
   cd $token
   echo $token
   echo $token.pdb > phenix.refine.txt
   phenix.ready_set $token.pdb
   FILE=$token.ligands.cif
   if [ -f "$FILE" ]
   then
     echo $token.updated.pdb *.mtz $token.ligands.cif >> phenix.refine.txt
     phenix.refine $token.updated.pdb *.mtz $token.ligands.cif
   else
     echo $token.pdb *.mtz >> phenix.refine.txt
     phenix.refine $token.pdb *.mtz
   fi
   mkdir WC
   mkdir  HG
   cp $token.update_refine_001.mtz WC
   cp $token.refine_001.mtz WC
   cp $token.update_refine_001.mtz HG
   cp $token.refine_001.mtz HG
   cd ..
 done
