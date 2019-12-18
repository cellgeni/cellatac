#!/bin/bash

set -euo pipefail

input=${1:?Need region bed file}
output=${2:?Need name for output file}

perl -ane 'local $"="_"; $i=$.-1; print "$i\t$F[0]:$F[1]-$F[2]\n"' $input > $output

#     perl -ane 'local $"="_"; print "$ENV{b}\t@F[0..2]\t$F[3]\n"' $f

