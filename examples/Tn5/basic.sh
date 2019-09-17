#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# standard input
echo; echo '#----- INPUT DNA: 5 ng ----#'
pyTecanFluent Tn5 --prefix /tmp/Tn5_basic_96well \
  $DIR/../../tests/data/basic_96well.txt

# lower input of DNA
echo; echo '#----- INPUT DNA: 1 ng -----#'
pyTecanFluent Tn5 --prefix /tmp/Tn5_basic_96well \
  --sample-conc 1 --sample-volume 1 \
  $DIR/../../tests/data/basic_96well.txt 

# higher input of DNA
echo; echo '#----- INPUT DNA: 10 ng -----#'
pyTecanFluent Tn5 --prefix /tmp/Tn5_basic_96well \
  --sample-conc 5 --sample-volume 2 \
  $DIR/../../tests/data/basic_96well.txt 
