#!/bin/bash

# This script reads many bedtools coverage output files and combines them in a matrix
# in mcl native binary format. This is very fast for reading, writing and subsetting.
# The most frequent windows are picked using this matrix, and a subsetted matrix
# in matrixmarket format is generated.

# TODO:
# make sure all barcodes have at least 1 window > 0.
# mcxdump is not guaranteed to output integers (for larger numbers) due to %g format (e.g. 9.84373e+06)
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

(cd $inputdir && ls -1 *.w5k.txt) | cut -f 1 -d '.' | nl -v0 -nln -w1 > bar.tab

export MCLXIOFORMAT=8   # force native binary format, it's 20-30 times faster.

if $force || ! mcx query --dim -imx bar2win.mcx; then

  for f in $inputdir/*.w5k.txt; do

    g=${f##*/}
    export b=${g%.w5k.txt}

    perl -ane 'local $"="_"; print "$ENV{b}\t@F[0..2]\t$F[3]\n"' $f

  done | \
  mcxload --stream-split -abc - -strict-tabc bar.tab -strict-tabr win.tab --write-binary -o bar2win.mcx
else
>&2 echo "Reusing bar2win.mcx"
fi

export MCLXIOVERBOSITY=2

>&2 echo "Computing transpose matrix"
  mcxi /bar2win.mcx lm tp /win2bar.mcx wm

>&2 echo "Computing barcode stats bar.stats"
  mcx query -imx bar2win.mcx -tab bar.tab > bar.stats
>&2 echo "Computing window stats win.stats, win$fTAG.stats, win$fTAG.names"
  mcx query -imx win2bar.mcx -tab win.tab > win.stats
  (sort -nrk 2 win.stats || true) | head -n $num_sites > win$fTAG.stats
  cut -f 1 win$fTAG.stats > win$fTAG.names

>&2 echo "Selecting windows using win$fTAG.names, win$fTAG.tab"
  perl -ane 'BEGIN{open(W,"<win$ENV{fTAG}.names")||die;%h=map{chomp;($_,1)}<W>}print if $h{$F[1]}' win.tab > win$fTAG.tab



>&2 echo "Computing submatrix using win$fTAG.tab, win2bar$fTAG.mcx"
  mcxsubs -imx win2bar.mcx --from-disk -tab win$fTAG.tab "dom(c, t()), out(win2bar$fTAG.mcx,wb)"

>&2 echo "Remapping indices, win2bar$fTAG.reindexed"
  mcxmap -imx win2bar$fTAG.mcx -make-mapc win2bar$fTAG.map -o win2bar$fTAG.reindexed

>&2 echo "Remapping tab file, win$fTAG.tab.reindexed"
  mcxmap -tab win$fTAG.tab -map win2bar$fTAG.map > win$fTAG.tab.reindexed

  cut -f 2 win$fTAG.tab.reindexed  > regions$num_sites.txt
  cut -f 2 bar.tab                > barcode.txt

  n_regions=$num_sites
  n_barcodes=$(cat barcode.txt | wc -l)
  n_entries=$(datamash sum 2 < win$fTAG.stats)

>&2 echo "Writing matrix market file"
(
cat <<EOH
%%MatrixMarket matrix coordinate integer general
$n_barcodes $n_regions $n_entries
EOH
mcxdump --transpose -imx win2bar$fTAG.reindexed | perl -ane 'BEGIN{$,="\t"} print $F[0]+1, $F[1]+1, "$F[2]\n";'
) | gzip > mtx.gz


