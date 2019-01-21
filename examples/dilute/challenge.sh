#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate & destination plates; >96 samples
INPUT=$DIR'/../../tests/data/dilutionfile_samples_Hagai.xlsx'
echo "Using input file: " $INPUT
pyTecanFluent dilute --prefix /tmp/TECAN_dilute_challenge $INPUT

