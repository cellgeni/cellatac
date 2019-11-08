#!/bin/bash

set -euo pipefail

inputfile=${1:?Need sam inputfile}

samtools view -b $inputfile -o ${inputfile%.sam}.bam

