#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate
pyTecanFluent LITE --prefix /tmp/LITE_1plate \
  --pcr-mm-labware-type Falcon_10ml_3X_tube_rack \
  $DIR/../../tests/data/LITE_2-input-plates.xlsx


