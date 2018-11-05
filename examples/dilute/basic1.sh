#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate & destination plates
INPUT=$DIR'/../../tests/data/conc_file1.txt'
echo "Using input file: " $INPUT
pyTecanFluent dilute --prefix /tmp/TECAN_dilute_conc1 $INPUT

# 96-well source plates & destination plates
INPUT=$DIR'/../../tests/data/conc_file2.txt'
echo "Using input file: " $INPUT
pyTecanFluent dilute --prefix /tmp/TECAN_dilute_conc2 $INPUT


# low concentrations
#INPUT=$DIR'/../../tests/data/conc_zach1.txt'
#echo "Using input file: " $INPUT
#pyTecanFluent dilute --prefix /tmp/TECAN_dilute_zach1 $INPUT
