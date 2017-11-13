#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate & destination plates; >96 samples
INPUT=$DIR'/../../tests/data/dilution_gt96.xlsx'
echo "Using input file: " $INPUT
pyTecanFluent dilute --prefix /tmp/TECAN_dilute_gt96 $INPUT

