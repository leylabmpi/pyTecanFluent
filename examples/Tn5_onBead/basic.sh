#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# standard input
echo && echo '#----- INPUT DNA: 5 ng ----#'
pyTecanFluent Tn5_onBead --prefix /tmp/Tn5-on-Bead_basic_96well \
  $DIR/../../tests/data/Tn5_onBead_96well.txt

