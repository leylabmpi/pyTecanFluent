#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate; dual barcoding
pyTecanFluent map2robot --prefix /ebio/abt3_projects/databases/TECAN/worklists/S5_N7 \
  $DIR/../../tests/data/samples_S5-N7.xlsx
