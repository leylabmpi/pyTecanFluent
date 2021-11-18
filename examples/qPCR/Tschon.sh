#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate & destination plates
INFILE=$DIR'/../../tests/data/Tschon_qPCR.xlsx'
echo "Using input file:" $INFILE
pyTecanFluent qPCR --prefix /tmp/TECAN_qPCR $INFILE

