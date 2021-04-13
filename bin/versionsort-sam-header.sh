#!/bin/bash

# Dual gex+atac sequencing delivers the contigs in lexicographic sort order.
# As we combine them with regular atac data and are used to version sort order,
# this script is used to re-sort the header of the sam file.  Note that we sort
# the fragment data as well in the get scripts (external to cellatac), using
#
# zcat $frag | sort --parallel=4 -k 1,1V -k 2,2n -k 3,3n | gzip > fragments.tsv.gz


set -euo pipefail

sam_in=${1?Please supply input sam (header) file name}
sam_out=${2?Please supply output sam (header) file name}

(
                                      # grab everything before the SQ section.
perl -pe 'exit if /^\@SQ/' $sam_in;
                                      # grab the SQ section and sort it.
grep "^@SQ" $sam_in | sort -V;
                                      # grab everything after the SQ section.
                                      # perl range operator .. is bistable.
perl -ne '$seen_sq = /^\@SQ/..1; $not_sq = !/^\@SQ/; print if $seen_sq && $not_sq' $sam_in
) > $sam_out

