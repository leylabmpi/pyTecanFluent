#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate
pyTecanFluent Tn5 --prefix /tmp/Tn5_T20 \
  $DIR/../../tests/data/Tn5_Taichi20.txt

