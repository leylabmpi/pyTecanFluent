#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate & destination plate; no mapping file
INPUT=$DIR'/../../tests/data/PCR-run1.xlsx'
echo "Using input file: " $INPUT
pyTecanFluent pool --prefix /tmp/TECAN_pool1 $INPUT

# 96-well source plate & destination plate; with mapping file
INPUT=$DIR'/../../tests/data/PCR-run1.xlsx'
echo "Using input file: " $INPUT
pyTecanFluent pool \
	      --prefix /tmp/TECAN_pool1 \
	      --mapfile $DIR'/../../tests/data/samples_S5-N7.xlsx' \
	      $INPUT
