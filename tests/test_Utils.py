#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import
## batteries
import os
import sys
import pytest
## 3rd party
import pandas as pd
## package
from pyTecanFluent import Utils


# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')

# tests
def test_make_range():
    x = Utils.make_range('all')
    assert x is None

    x = Utils.make_range('0')
    assert x == [0]
    
    x = Utils.make_range('1,2,5')
    assert x == [1,2,5]
    
    x = Utils.make_range('1,2,5-6')
    assert x == [1,2,5,6]

    x = Utils.make_range('1-3,5-6')
    assert x == [1,2,3,5,6]
    
def test_range_zeroindex():
    x = Utils.make_range('all', set_zero_index=True)
    assert x is None

    with pytest.raises(ValueError):
        Utils.make_range('0', set_zero_index=True)

    x = Utils.make_range('1,2,5', set_zero_index=True)
    assert x == [0,1,4]

    x = Utils.make_range('1,2,5-6', set_zero_index=True)
    assert x == [0,1,4,5]
        
def test_check_gwl():
    gwl_file = os.path.join(data_dir, 'multi_dispense.gwl')
    ret = Utils.check_gwl(gwl_file)
    assert ret is None

