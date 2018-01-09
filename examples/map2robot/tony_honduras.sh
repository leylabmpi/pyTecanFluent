#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well destination plate, subset of total samples
pyTecanFluent map2robot --prefix /tmp/tony_honduras_small \
  $DIR/../../tests/data/tony_honduras_map_small.tsv

# 96-well destination plate
pyTecanFluent map2robot --prefix /tmp/tony_honduras \
  $DIR/../../tests/data/tony_honduras_map_plate1.tsv
