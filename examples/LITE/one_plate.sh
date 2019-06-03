#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate
pyTecanFluent LITE --prefix /tmp/LITE_1plate \
	      --pcr-mm-volume 7 \
	      --pcr-mm-labware-type "10ml Falcon" \
	      $DIR/../../tests/data/LITE_1plate.xlsx


# # 384-well source plate
# pyTecanFluent map2robot --prefix /tmp/Map2Robot_basic_384well \
#   $DIR/../../tests/data/basic_384well.txt

# # 96-well source plate, 384-well destination plate
# pyTecanFluent map2robot --prefix /tmp/Map2Robot_basic_96-384well \
#   --dest-type '384 Well Biorad PCR' \
#   $DIR/../../tests/data/basic_96well.txt

# # 96-well source plate, 96-well designation, no primers
# pyTecanFluent map2robot --prefix /tmp/Map2Robot_basic_96well_noPrimer \
#   --prm-volume 0 \
#   $DIR/../../tests/data/basic_96well.txt
