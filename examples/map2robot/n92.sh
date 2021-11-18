#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate
pyTecanFluent map2robot --prefix /tmp/Map2Robot_n92 \
  $DIR/../../tests/data/map2robot_n92.tsv

