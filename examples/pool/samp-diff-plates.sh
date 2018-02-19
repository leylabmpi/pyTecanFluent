#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate & destination plate; no mapping file; source samples on different plates
INPUT=$DIR'/../../tests/data/TECAN_pool_samples-diff-plates.txt'
echo "Using input file: " $INPUT
pyTecanFluent pool --prefix /tmp/TECAN_pool1 $INPUT

