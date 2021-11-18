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
from pyTecanFluent import Fluent
from pyTecanFluent import Labware

# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')

# tests
def test_db():
    db = Fluent.db()    
    RackTypes = db.RackTypes()
    assert isinstance(RackTypes, list)
    
    RackType = db.RackTypes()[0]
    v = db.get_labware(RackType)
    assert isinstance(v, dict)
    
def test_aspirate():
    asp = Fluent.Aspirate()
    asp.RackLabel = 'test'
    asp.RackType = 'test'
    assert isinstance(asp.cmd(), str)

def test_dispense():
    disp = Fluent.Dispense()
    disp.RackLabel = 'test'
    disp.RackType = 'test'
    assert isinstance(disp.cmd(), str)

def test_comment():
    c = Fluent.Comment()
    assert isinstance(c.cmd(), str)

def test_waste():
    w = Fluent.Waste()
    assert isinstance(w.cmd(), str)
    
def test_gwl():
    gwl = Fluent.gwl()
    with pytest.raises(AssertionError):
        gwl.add(Fluent.Aspirate())

def test_labware():
    lw = Labware.labware()
    gwl = Fluent.gwl()
    ret = lw.add_gwl(gwl)
    assert ret is None
        
