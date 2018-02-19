#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate & destination plate; no mapping file; dest well offset
INPUT=$DIR'/../../tests/data/Pooling_well_offset.xlsx'
echo "Using input file: " $INPUT
pyTecanFluent pool --dest-start 71 --prefix /tmp/TECAN_pool1 $INPUT 
