#!/bin/bash

# This script reads many bedtools coverage peak output files and outputs
# them in matrixmarket format. See ca_top_regions.sh for more comments.

set -euo pipefail

peakbedfile=
filelist=
cellnames=
force=false


while getopts :c:i:p:Fh opt
do
    case "$opt" in
    c)
      cellnames=$OPTARG
      ;;
    i)
      filelist=$OPTARG
      ;;
    p)
      peakbedfile=$OPTARG
      ;;
    F)
      force=true
      ;;
    h)
      cat <<EOU
-c cell ID file (barcodes usually; one id per line)
-w window file
-i file with file locations inside, one file per line, full path name
EOU
      exit
      ;;
    :) echo "Flag $OPTARG needs argument"
        exit 1;;
    ?) echo "Flag $OPTARG unknown"              # Yes that's right, $OPTARG. bash bish bosh.
        exit 1;;
   esac
done

if [[ -z $filelist || -z $peakbedfile ]]; then
   echo "Need -i filelistfile and -p peak bed file! (see -h)"
   false
fi


ca_region_maketab.sh $peakbedfile __peak.tab

nl -v0 -nln -w1 < $cellnames > __cell.tab

cut -f 2 __peak.tab > peaks.txt
cut -f 2 __cell.tab > cells.txt         # will be identical to $cellnames.


export MCLXIOFORMAT=8   # force native binary format, it's 20-30 times faster.

# # ###################################
 #  Stream all files for matrix loading
#   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
if $force || [[ ! -e peak2cell.mcx ]]; then

  while read f; do

    g=${f##*/}
    export b=${g%.mp.txt}

    perl -ane 'local $"="_"; print "$F[0]:$F[1]-$F[2]\t$ENV{b}\t$F[3]\n"' $f

  done < "$filelist" | mcxload \
          --stream-split -abc - -strict-tabr __cell.tab -strict-tabc __peak.tab --write-binary -o peak2cell.mcx
else
>&2 echo "Reusing peak2cell.mcx"
fi


mcx query -imx peak2cell.mcx -o peak2cell.stats -tab __peak.tab
n_entries=$(tail -n +2 peak2cell.stats | perl -ane '$S+=$F[1]; END{print "$S\n";}')


# # #############################
 #  Output to matrixmarket format
#   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ca_make_mmtx.sh -r peaks.txt -c cells.txt -m peak2cell.mcx -e $n_entries -t integer -o peaks_bc_matrix.mmtx.gz
ca_make_mmtx.sh -r peaks.txt -c cells.txt -m peak2cell.mcx -e $n_entries -t integer -o bc_peaks_matrix.mmtx.gz -T
# various naming schemes exist. Maybe we'll use bc.
ln cells.txt bc.txt


