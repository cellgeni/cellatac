#!/bin/bash

set -euo pipefail

FILE=$1

sort -k1,1V -k2,2n -k 3,3n $FILE > t
rm $FILE
mv t $FILE
