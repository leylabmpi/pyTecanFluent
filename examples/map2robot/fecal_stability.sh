#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well destination plate
pyTecanFluent map2robot --prefix /tmp/fecal_stab_96well \
  $DIR/../../tests/data/mapping_file_fecal_stability.txt

# 394-well destination plate
pyTecanFluent map2robot --prefix /tmp/fecal_stab_384well \
	      --dest-type '384 Well Biorad PCR' \
	      --mm-one-source \
	      $DIR/../../tests/data/mapping_file_fecal_stability.txt
