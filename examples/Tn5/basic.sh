#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate
pyTecanFluent Tn5 --prefix /tmp/Tn5_basic_96well \
  $DIR/../../tests/data/basic_96well.txt

