#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import
## batteries
import os
import sys
import shutil
import pytest
## 3rd party
import pandas as pd
## package
from pyTecanFluent import Tn5
from pyTecanFluent import Utils


# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
def test_help(script_runner):
    ret = script_runner.run('pyTecanFluent', 'Tn5', '-h')
    assert ret.success

def test_basic(script_runner, tmp_path):
    output_prefix = os.path.join(str(tmp_path), 'basic')
    conc_file = os.path.join(data_dir, 'basic_96well.txt')
    ret = script_runner.run('pyTecanFluent', 'Tn5', '--prefix',
                            output_prefix, conc_file)
    assert ret.success
    
def test_T20(script_runner, tmp_path):
    output_prefix = os.path.join(str(tmp_path), 'Taichi20')
    conc_file = os.path.join(data_dir, 'Tn5_Taichi20.txt')
    ret = script_runner.run('pyTecanFluent', 'Tn5', '--prefix',
                            output_prefix, conc_file)
    assert ret.success
    
def test_T60(script_runner, tmp_path):
    output_prefix = os.path.join(str(tmp_path), 'Taichi60')
    conc_file = os.path.join(data_dir, 'Tn5_Taichi60.txt')
    ret = script_runner.run('pyTecanFluent', 'Tn5',
                            '--tag-n-tip-reuse', '4',
                            '--prefix', output_prefix, conc_file)
    assert ret.success
    
