#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate & destination plates
INFILE=$DIR'/../../tests/data/2019_06_02_CM12_Flynn_384.xlsx'
echo "Using input file:" $INFILE
pyTecanFluent qPCR --mm-type "10ml Falcon" --prefix /tmp/TECAN_qPCR $INFILE

