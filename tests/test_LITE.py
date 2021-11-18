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
from pyTecanFluent import LITE
from pyTecanFluent import Utils


# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
def test_help(script_runner):
    ret = script_runner.run('pyTecanFluent', 'LITE', '-h')
    assert ret.success

def test_basic(script_runner, tmp_path):
    output_prefix = os.path.join(str(tmp_path), 'basic')
    conc_file = os.path.join(data_dir, 'basic_96well.txt')
    ret = script_runner.run('pyTecanFluent', 'LITE', '--prefix',
                            output_prefix, conc_file)
    assert ret.success
    
def test_onePlate(script_runner, tmp_path):
    output_prefix = os.path.join(str(tmp_path), 'onePlate')
    conc_file = os.path.join(data_dir, 'LITE_1plate.xlsx')
    ret = script_runner.run('pyTecanFluent', 'LITE',
                            '--pcr-mm-volume', '7',
                            '--pcr-mm-labware-type', '10ml Falcon',
                            '--prefix', output_prefix, conc_file)
    assert ret.success

