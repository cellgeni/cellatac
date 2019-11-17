#!/bin/bash

# This script reads many bedtools coverage peak output files and outputs
# them in matrixmarket format. See cellatac_top_regions.sh for more comments.

set -euo pipefail

windowfile=
filelist=
cellnames=
force=false


while getopts :c:i:w:Fh opt
do
    case "$opt" in
    c)
      cellnames=$OPTARG
      ;;
    i)
      filelist=$OPTARG
      ;;
    w)
      windowfile=$OPTARG
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

if [[ -z $filelist || -z $windowfile ]]; then
   echo "Need -i filelistfile and -w window file! (see -h)"
   false
fi


cellatac_region_maketab.sh $windowfile win.tab

nl -v0 -nln -w1 < $cellnames > cell.tab

cut -f 2 win.tab > peaks.txt
cut -f 2 cell.tab > cells.txt         # will be identical to $cellnames.


export MCLXIOFORMAT=8   # force native binary format, it's 20-30 times faster.

# # ###################################
 #  Stream all files for matrix loading
#   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
if $force || ! -e peak2cell.mcx; then

  while read f; do

    g=${f##*/}
    export b=${g%.mp.txt}

    perl -ane 'local $"="_"; print "@F[0..2]\t$ENV{b}\t$F[3]\n"' $f

  done < "$filelist" | mcxload \
          --stream-split -abc - -strict-tabr cell.tab -strict-tabc win.tab --write-binary -o peak2cell.mcx
else
>&2 echo "Reusing peak2cell.mcx"
fi


mcx query -imx peak2cell.mcx -o peak2cell.stats
n_entries=$(tail -n +2 peak2cell.stats | perl -ane '$S+=$F[1]; END{print "$S\n";}')


# # #############################
 #  Output to matrixmarket format
#   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

cellatac_make_mmtx.sh -r peaks.txt -c cells.txt -m peak2cell.mcx -e $n_entries -t integer -o cell2peak.gz


