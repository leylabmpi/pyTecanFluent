#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate & destination plates
#INFILE=$DIR'/../../tests/data/qPCR_setup/qPCR_Zach_plate1.xlsx'
INFILE=$DIR'/../../tests/data/Cmin_qPCR.tsv'
echo "Using input file:" $INFILE
pyTecanFluent qPCR \
  --prefix /tmp/TECAN_qPCR $INFILE \
  --mm-type "5ml Eppendorf waste"

