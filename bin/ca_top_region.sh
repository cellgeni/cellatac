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

wintab=
force=false
num_sites=20000
celltab=
cell2winmtx=

o_winstats=win.stats
o_cellstats=cell.stats
o_regions_names=regions.names
o_cells_names=cells.names


while getopts :m:w:N:W:c:C:R:n:Fh opt
do
    case "$opt" in
    c)
      celltab=$OPTARG
      ;;
    C)
      o_cellstats=$OPTARG
      ;;
    w)
      wintab=$OPTARG
      ;;
    W)
      o_winstats=$OPTARG
      ;;
    N)
      o_cells_names=$OPTARG
      ;;
    R)
      o_regions_names=$OPTARG
      ;;
    m)
      cell2winmtx=$OPTARG
      ;;
    n)
      num_sites=$OPTARG
      ;;
    F)
      force=true
      ;;
    h)
      cat <<EOU
-c cell tab file
-w window tab file
-C cell stats (output) file
-W window stats (output) file
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

if [[ -z $cell2winmtx || -z $celltab || -z $wintab || -z $o_winstats || -z $o_cellstats ]]; then
   echo "Need -m matrixfile -c celltab and -w wintab -C cellstats -W winstats ! (see -h)"
   false
fi


export MCLXIOFORMAT=8   # force native binary format, it's 20-30 times faster.


# # ###############################################################
 #  Compute cell/region statistics both ways and select top regions
#   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
export MCLXIOVERBOSITY=2

>&2 echo "Computing transpose matrix"
  mcxi /$cell2winmtx lm tp /__win2cell.mcx wm

>&2 echo "Computing cell stats __cell.stats"
  mcx query -imx $cell2winmtx -tab $celltab > __cell.stats
>&2 echo "Computing window stats $o_winstats, __win_TAG.stats, __win_TAG.names"
  mcx query -imx __win2cell.mcx -tab $wintab > $o_winstats
  pivot=$((sort -nrk 2 $o_winstats || true) | head -n $num_sites | tail -n 1 | cut -f 2)
  perl -ane 'print if $F[1] >= '$pivot $o_winstats > __win_TAG.stats
  cut -f 1 __win_TAG.stats > __win_TAG.names

>&2 echo "Selecting windows using __win_TAG.names, __win_TAG.tab"
  perl -ane 'BEGIN{open(W,"<__win_TAG.names")||die;%h=map{chomp;($_,1)}<W>}print if $h{$F[1]}' $wintab > __win_TAG.tab


# # ################################
 #  Compute submatrix, remap indices
#   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

>&2 echo "Computing submatrix using __win_TAG.tab, __win2cell_TAG.mcx"
  mcxsubs -imx __win2cell.mcx --from-disk -tab __win_TAG.tab "dom(c, t()), out(__win2cell_TAG.mcx,wb)"

>&2 echo "Remapping indices, __win2cell_TAG.reindexed"
  mcxmap -imx __win2cell_TAG.mcx -make-mapc __win2cell_TAG.map -o __win2cell_TAG.reindexed

>&2 echo "Remapping tab file, __win_TAG.tab.reindexed"
  mcxmap -tab __win_TAG.tab -map __win2cell_TAG.map > __win_TAG.tab.reindexed

>&2 echo "Transposing matrix, __cell2win_TAG.mcx"
  mcxi /__win2cell_TAG.reindexed lm tp /__cell2win_TAG.mcx wm

>&2 echo "Computing cell stats files on reduced region set, __cell_TAG.stats"
  mcx query -imx __cell2win_TAG.mcx -tab $celltab > __cell_TAG.stats
  ln __cell_TAG.stats $o_cellstats

  cut -f 2 __win_TAG.tab.reindexed  > $o_regions_names
  cut -f 2 $celltab                 > $o_cells_names


# # #############################
 #  Output to matrixmarket format
#   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

>&2 echo "Producing matrixmarket format"
  n_entries=$(perl -ane '$S+=$F[1]; END{print "$S\n";}' __win_TAG.stats)

  ca_make_mmtx.sh -r __win_TAG.names -c $o_cells_names -m __win2cell_TAG.reindexed -e $n_entries -t pattern -o mtx.gz


