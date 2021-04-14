#!/bin/bash

set -euo pipefail

manifest=${1?Please provide manifest file}

totalcounts=0
while read sampleid path;
  do
    thefile="$path/singlecell.csv"
    if [[ ! -e $thefile ]]
    then 
      echo "singlecell.csv not present for sample" $sampleid && false
    else
      counts=$(get-col.py -i $thefile -s "," -c "is__cell_barcode" -f "1" -o "," | wc -l)
      echo $sampleid $counts
      ((totalcounts+=counts))
    fi;
done < <(cut -f 2,3 $manifest)
echo "total:" $totalcounts
