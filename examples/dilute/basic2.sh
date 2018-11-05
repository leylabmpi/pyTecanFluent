#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate & destination plates
INPUT=$DIR'/../../tests/data/conc_hagay1.txt'
echo "Using input file: " $INPUT
pyTecanFluent dilute --min-volume 2 --max-volume 20 --min-total 20 --dilution 0.5 --prefix /tmp/TECAN_dilute_conc1 $INPUT


