#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import
## batteries
import os
import sys
import shutil
import tempfile
import pytest
#import unittest
## 3rd party
import pandas as pd
## package
from pyTecanFluent import Tn5_onBead
from pyTecanFluent import Utils


# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
def test_help(script_runner):
    ret = script_runner.run('pyTecanFluent', 'Tn5_onBead', '-h')
    assert ret.success

def test_basic(script_runner, tmp_path):
    output_prefix = os.path.join(str(tmp_path), 'basic')
    conc_file = os.path.join(data_dir, 'Tn5_onBead_96well.txt')
    ret = script_runner.run('pyTecanFluent', 'Tn5_onBead', '--prefix',
                            output_prefix, conc_file)
    assert ret.success

    
