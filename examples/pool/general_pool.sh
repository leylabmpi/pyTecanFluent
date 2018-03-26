#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 384+96-well source plate & destination plate; no mapping file; dest = 1 tube
INPUT=$DIR'/../../tests/data/pooling_general.xlsx'
echo "Using input file: " $INPUT
pyTecanFluent pool --prefix /tmp/TECAN_pool_general $INPUT --dest-type "2ml Eppendorf"


# 384+96-well source plate & destination plate; no mapping file; dest = 1 tube; differing volumes
INPUT=$DIR'/../../tests/data/pooling_general_diff-volumes.xlsx'
echo "Using input file: " $INPUT
pyTecanFluent pool --prefix /tmp/TECAN_pool_general $INPUT --dest-type "2ml Eppendorf" --volume-col "Volume"
