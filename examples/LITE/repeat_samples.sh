#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 96-well source plate
pyTecanFluent LITE --prefix /tmp/LITE_repeat \
	      --tag-mm-volume 0 --pcr-mm-volume 0 \
	      --dest-type "PCR Adapter 96 Well and 96 Well Eppendorf TwinTec PCR" \
	      $DIR/../../tests/data/Gewirtz_samples_repeat_LITE.xlsx


