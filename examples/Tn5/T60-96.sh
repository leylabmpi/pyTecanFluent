#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 60 samples
pyTecanFluent Tn5 --prefix /tmp/Tn5_T60 \
      --tag-n-tip-reuse 4 \
      $DIR/../../tests/data/Tn5_Taichi60.txt

# 96 samples
pyTecanFluent Tn5 --prefix /tmp/Tn5_T96 \
      --tag-n-tip-reuse 6 \
      $DIR/../../tests/data/Tn5_Taichi96.txt
