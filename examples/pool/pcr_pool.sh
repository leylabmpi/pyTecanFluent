#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate & destination plate; no mapping file
INPUT=$DIR'/../../tests/data/pooling_16Spcr.xls'
echo "Using input file: " $INPUT
pyTecanFluent pool --prefix /tmp/TECAN_pool_pcr $INPUT

