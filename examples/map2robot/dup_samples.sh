#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate
pyTecanFluent map2robot --prefix /tmp/Map2Robot_basic_dupSamples \
  $DIR/../../tests/data/TECAN_Georg_PCR_dupSamples.xlsx


