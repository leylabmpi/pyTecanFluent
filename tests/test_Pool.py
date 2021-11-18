#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import
## batteries
import os
import sys
import shutil
import tempfile
import pytest
## 3rd party
import pandas as pd
## package
from pyTecanFluent import Pool


# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
def test_noMap(script_runner, tmp_path):    
    output_prefix = os.path.join(str(tmp_path), 'noMap')
    pcr_file = os.path.join(data_dir, 'PCR-run1.xlsx')
    ret = script_runner.run('pyTecanFluent', 'pool',                            
                            '--prefix', output_prefix, pcr_file)
    assert ret.success

def test_map(script_runner, tmp_path):
    output_prefix = os.path.join(str(tmp_path), 'map')
    pcr1_file = os.path.join(data_dir, 'PCR-run1.xlsx')
    pcr2_file = os.path.join(data_dir, 'PCR-run2.xlsx')
    map_file = os.path.join(data_dir, 'samples_S5-N7.xlsx')
    ret = script_runner.run('pyTecanFluent', 'pool',
                            '--prefix', output_prefix,
                            '--mapfile', map_file,
                            pcr1_file, pcr2_file)
    assert ret.success
    
