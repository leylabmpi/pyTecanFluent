#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate
pyTecanFluent map2robot --prefix /tmp/Map2Robot_basic_96well \
  $DIR/../../tests/data/basic_96well.txt

# 384-well source plate
pyTecanFluent map2robot --prefix /tmp/Map2Robot_basic_384well \
  $DIR/../../tests/data/basic_384well.txt

# 96-well source plate, 384-well destination plate
pyTecanFluent map2robot --prefix /tmp/Map2Robot_basic_96-384well \
  --dest-type '384 Well Biorad PCR' \
  $DIR/../../tests/data/basic_96well.txt
