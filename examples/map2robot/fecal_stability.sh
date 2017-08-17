#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well destination plate
TECAN map2robot --prefix /ebio/abt3_projects/databases/TECAN/worklists/fecal_stab_96well \
  $DIR/../../tests/data/mapping_file_fecal_stability.txt

# 394-well destination plate
TECAN map2robot --prefix /ebio/abt3_projects/databases/TECAN/worklists/fecal_stab_384well \
  --desttype 384 \
  $DIR/../../tests/data/mapping_file_fecal_stability.txt
