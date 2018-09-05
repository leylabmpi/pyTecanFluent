#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate, 96-well designation, no primers
pyTecanFluent map2robot --prefix /tmp/Map2Robot_basic_96well_noPrimer \
  --prm-volume 0 \
  $DIR/../../tests/data/basic_96well.txt

# 96-well source plate, 96-well designation, no primers, no water
pyTecanFluent map2robot --prefix /tmp/Map2Robot_basic_96well_noPrimer-water \
  --prm-in-mm --water-in-mm \
  $DIR/../../tests/data/basic_96well.txt

# 96-well source plate, 96-well designation, no mastermix reagent distribution
pyTecanFluent map2robot --prefix /tmp/Map2Robot_basic_96well_noMM-RD \
  --n-multi-disp 1 --mm-liq "MasterMix Free Single Wall Disp" \
  $DIR/../../tests/data/basic_96well.txt
