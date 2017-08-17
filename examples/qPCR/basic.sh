#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate & destination plates
INFILE=$DIR'/../../tests/data/qPCR_setup/qPCR_Zach_plate1.xlsx'
echo "Using input file:" $INFILE
TECAN qPCR --prefix /tmp/TECAN_qPCR $INFILE

