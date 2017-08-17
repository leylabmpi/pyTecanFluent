#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate & destination plates
INPUT=$DIR'/../../tests/data/conc_file1.txt'
echo "Using input file: " $INPUT
TECAN dilute --prefix /tmp/TECAN_dilute $INPUT

