#!/bin/bash

# This script reads many bedtools coverage output files and combines them in a
# matrix in mcl native binary format. This format is very fast for reading,
# writing and subsetting.  The most frequent windows are picked and a subsetted
# matrix in matrixmarket format is generated.

# TODO:
# make sure all barcodes have at least 1 window > 0.
# mcxdump is not guaranteed to output integers (for larger numbers) due to %g format (e.g. 9.84373e+06)
# -> not relevant really, as currently approach normalises everything to 1.
# file existence test currently present can cause resumption to fail if corrupted (remove or add force)

set -euo pipefail

filelist=
windowfile=
force=false
num_sites=20000
cellnames=
fTAG=


while getopts :w:c:i:I:n:Fh opt
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
    n)
      num_sites=$OPTARG
      ;;
    F)
      force=true
      ;;
    h)
      cat <<EOU
-c cell ID file (barcodes usually; one id per line)
-w window file
-i file with file locations inside, one file per line, full path name
-F force first time-consuming step even if result file is present
-n number of sites to consider (default $num_sites)
EOU
      exit
      ;;
    :) echo "Flag $OPTARG needs argument"
        exit 1;;
    ?) echo "Flag $OPTARG unknown"              # Yes that's right, $OPTARG. bash bish bosh.
        exit 1;;
   esac
done

if [[ -z $filelist || -z $windowfile= ]]; then
   echo "Need -i filelistfile and -w windowfile=! (see -h)"
   false
fi

if [[ -z $fTAG ]]; then
  fTAG=$num_sites
fi
export fTAG

ca_region_maketab.sh $windowfile win.tab
# perl -ane 'local $"="_"; $i=$.-1; print "$i\t@F[0..2]\n"' $windowfile > win.tab

# (cd $inputdir && ls -1 *.w5k.txt) | cut -f 1 -d '.' | nl -v0 -nln -w1 > cell.tab
nl -v0 -nln -w1 < $cellnames > cell.tab

export MCLXIOFORMAT=8   # force native binary format, it's 20-30 times faster.

# # ###################################
 #  Stream all files for matrix loading
#   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
if $force || ! -e cell2win.mcx; then

  # for f in $inputdir/*.w5k.txt; do
  while read f; do

    g=${f##*/}
    export b=${g%.w5k.txt}

    perl -ane 'local $"="_"; print "$ENV{b}\t@F[0..2]\t$F[3]\n"' $f

  done < "$filelist" | mcxload \
          --stream-split -abc - -strict-tabc cell.tab -strict-tabr win.tab --write-binary -o cell2win.mcx
else
>&2 echo "Reusing cell2win.mcx"
fi


# # ###############################################################
 #  Compute cell/region statistics both ways and select top regions
#   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
export MCLXIOVERBOSITY=2

>&2 echo "Computing transpose matrix"
  mcxi /cell2win.mcx lm tp /win2cell.mcx wm

>&2 echo "Computing cell stats cell.stats"
  mcx query -imx cell2win.mcx -tab cell.tab > cell.stats
>&2 echo "Computing window stats win.stats, win$fTAG.stats, win$fTAG.names"
  mcx query -imx win2cell.mcx -tab win.tab > win.stats
  pivot=$((sort -nrk 2 win.stats || true) | head -n $num_sites | tail -n 1 | cut -f 2)
  perl -ane 'print if $F[1] >= '$pivot win.stats > win$fTAG.stats
  cut -f 1 win$fTAG.stats > win$fTAG.names

>&2 echo "Selecting windows using win$fTAG.names, win$fTAG.tab"
  perl -ane 'BEGIN{open(W,"<win$ENV{fTAG}.names")||die;%h=map{chomp;($_,1)}<W>}print if $h{$F[1]}' win.tab > win$fTAG.tab


# # ################################
 #  Compute submatrix, remap indices
#   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

>&2 echo "Computing submatrix using win$fTAG.tab, win2cell$fTAG.mcx"
  mcxsubs -imx win2cell.mcx --from-disk -tab win$fTAG.tab "dom(c, t()), out(win2cell$fTAG.mcx,wb)"

>&2 echo "Remapping indices, win2cell$fTAG.reindexed"
  mcxmap -imx win2cell$fTAG.mcx -make-mapc win2cell$fTAG.map -o win2cell$fTAG.reindexed

>&2 echo "Remapping tab file, win$fTAG.tab.reindexed"
  mcxmap -tab win$fTAG.tab -map win2cell$fTAG.map > win$fTAG.tab.reindexed

>&2 echo "Transposing matrix, cell2win$fTAG.mcx"
  mcxi /win2cell$fTAG.reindexed lm tp /cell2win$fTAG.mcx wm

>&2 echo "Computing cell stats files on reduced region set, cell$fTAG.stats"
  mcx query -imx cell2win$fTAG.mcx -tab cell.tab > cell$fTAG.stats

  cut -f 2 win$fTAG.tab.reindexed  > regions$num_sites.txt
  cut -f 2 cell.tab                > cells.txt


# # #############################
 #  Output to matrixmarket format
#   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

>&2 echo "Producing matrixmarket format"
  n_entries=$(perl -ane '$S+=$F[1]; END{print "$S\n";}' win$fTAG.stats)

  ca_make_mmtx.sh -r win$fTAG.names -c cells.txt -m win2cell$fTAG.reindexed -e $n_entries -t pattern -o mtx.gz

