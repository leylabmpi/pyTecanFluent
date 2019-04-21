#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate & destination plates; >96 samples
INPUT=$DIR'/../../tests/data/dilution_3plates.xlsx'
echo "Using input file: " $INPUT
pyTecanFluent dilute --prefix /tmp/TECAN_dilute_plate-order --dilution 1 $INPUT

