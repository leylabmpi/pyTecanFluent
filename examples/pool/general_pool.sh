#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate & destination plate; no mapping file; dest well offset
INPUT=$DIR'/../../tests/data/pooling_general.xlsx'
echo "Using input file: " $INPUT
pyTecanFluent pool --prefix /tmp/TECAN_pool_general $INPUT --dest-type "2ml Eppendorf"
