#!/bin/bash

set -euo pipefail

samp=${1:?Need sample ID}
export barc=${2:?Need barcode}

dataroot=/fix/this/script
dest=/fix/it/some/more

samtools view -b -f 2 -F 512 -F 256 $dataroot/$samp.cram | \
  bedtools bamtobed       | \
  cut -f 1,2,3            | \
  sort -k 1,1V -k 2,2n -k 3,3n    | \
  grep -i 'chr[a-z0-9][a-z0-9]*\>' | \
  uniq                    | \
  perl -ane 'chomp; print "$_\t$ENV{barc}\t1\n"' | \
  gzip                      \
  > $dest/$samp.frags.gz

