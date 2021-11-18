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
from pyTecanFluent import Fluent
from pyTecanFluent import Map2Robot
from pyTecanFluent import Utils


# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
def test_help(script_runner):
    ret = script_runner.run('pyTecanFluent', 'map2robot', '-h')
    assert ret.success

def test_basic(script_runner, tmp_path):
    output_prefix = os.path.join(str(tmp_path), 'basic')
    map_file = os.path.join(data_dir, 'basic_96well.txt')
    ret = script_runner.run('pyTecanFluent', 'map2robot', '--prefix',
                            output_prefix, map_file)
    assert ret.success

def test_single_barcode(script_runner, tmp_path):
    output_prefix = os.path.join(str(tmp_path), 'single-barcode')
    map_file = os.path.join(data_dir, 'mapping_file_fecal_stability.txt')
    ret = script_runner.run('pyTecanFluent', 'map2robot', '--prefix',
                            output_prefix, map_file)
    assert ret.success
        
