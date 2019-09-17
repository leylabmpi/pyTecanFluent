#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# standard input
pyTecanFluent Tn5 --prefix /tmp/Tn5_basic_96well \
  $DIR/../../tests/data/basic_96well.txt

# lower input of DNA
pyTecanFluent Tn5 --prefix /tmp/Tn5_basic_96well \
  --sample-conc 1 --sample-volume 1 \
  $DIR/../../tests/data/basic_96well.txt 
