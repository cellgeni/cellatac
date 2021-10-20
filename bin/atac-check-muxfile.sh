#!/bin/bash

muxfile=${1:?Need file path to manifest file, e.g. /some/path/manifest.txt}
output=${2:-chrom_file.txt}

> $output

cat $muxfile | while read i j k; do
  bcfile=$k/fragments.tsv.gz
  zcat $bcfile | cut -f 1 | uniq | grep -v "#" | while read l; do
    if [[ $l =~ 'chr' ]]; then
      echo "$j	$l" >> $output
    else
      echo "$j	$l" >> $output
    fi
    done;
done

if cut -f 2 $output | grep -q -w "^chr[0-9A-Z][0-9A-Z]*" && cut -f 2 $output | grep -q -w "^[0-9A-Z][0-9A-Z]*"; then
  echo "chromosome numbering is not consistent"
else
  echo "chromosome numbering is consistent"
fi
