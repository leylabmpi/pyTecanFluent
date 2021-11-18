#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well destination plate, subset of total samples
pyTecanFluent map2robot --prefix /tmp/PC_index \
  $DIR/../../tests/data/Map2Robot_Plate1_MITAB_PCR1_64_samples.xlsx

