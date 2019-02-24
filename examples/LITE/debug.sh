#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate
pyTecanFluent LITE --prefix /tmp/LITE_debug \
	      --tag-mm-volume 0 \
	      --sample-volume 0 \
	      $DIR/../../tests/data/silke_LITE_noSample.txt \

