#!/bin/bash

# This script reads many bedtools coverage output files and combines them in a
# matrix in mcl native binary format. This format very fast for reading,
# writing and subsetting.  The most frequent windows are picked and a subsetted
# matrix in matrixmarket format is generated.

# TODO:
# make sure all barcodes have at least 1 window > 0.
# mcxdump is not guaranteed to output integers (for larger numbers) due to %g format (e.g. 9.84373e+06)
# -> not relevant really, as currently approach normalises everything to 1.
# file existence test currently presents can cause resumption to fail if corrupted (remove or add force)
# -> added mcx query test, but still not sufficient.
# barcode index creation depends on '.' field separation in file name.

set -euo pipefail

inputdir=
genomefile=
force=false
num_sites=20000
fTAG=


while getopts :g:i:n:Fh opt
do
    case "$opt" in
    i)
      inputdir=$OPTARG
      ;;
    g)
      genomefile=$OPTARG
      ;;
    n)
      num_sites=$OPTARG
      ;;
    F)
      force=true
      ;;
    h)
      cat <<EOU
-g genome window file
-i input directory with bed coverage files
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

if [[ -z $inputdir || -z $genomefile ]]; then
   echo "Need -i inputdir and -g genomefile! (see -h)"
   false
fi

if [[ -z $fTAG ]]; then
  fTAG=$num_sites
fi
export fTAG

perl -ane 'local $"="_"; $i=$.-1; print "$i\t@F[0..2]\n"' $genomefile > win.tab

(cd $inputdir && ls -1 *.w5k.txt) | cut -f 1 -d '.' | nl -v0 -nln -w1 > cell.tab

export MCLXIOFORMAT=8   # force native binary format, it's 20-30 times faster.

# # ###################################
 #  Stream all files for matrix loading
#   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
if $force || ! -e cell2win.mcx; then

  for f in $inputdir/*.w5k.txt; do

    g=${f##*/}
    export b=${g%.w5k.txt}

    perl -ane 'local $"="_"; print "$ENV{b}\t@F[0..2]\t$F[3]\n"' $f

  done | \
  mcxload --stream-split -abc - -strict-tabc cell.tab -strict-tabr win.tab --write-binary -o cell2win.mcx
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
  n_regions=$(cat win$fTAG.names | wc -l)
  n_cells=$(cat cells.txt | wc -l)
  n_entries=$(perl -ane '$S+=$F[1]; END{print "$S\n";}' win$fTAG.stats)

>&2 echo "Writing matrix market file"
(
cat <<EOH
%%MatrixMarket matrix coordinate pattern general
$n_regions $n_cells $n_entries
EOH
mcxdump -imx win2cell$fTAG.reindexed --no-values | perl -ane 'print $F[0]+1, "\t", $F[1]+1, "\n";'
) | gzip > mtx.gz

