#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import
## batteries
import os
import sys
import shutil
import tempfile
import unittest
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
class Test_Map2Robot_main_basic(unittest.TestCase):

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.prefix = os.path.join(self.tmp_dir, 'output_')
        mapfile = os.path.join(data_dir, 'basic_96well.txt')
        self.args = Map2Robot.parse_args(['--prefix', self.prefix, mapfile])
        self.files = Map2Robot.main(self.args)

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)        

    def test_main_gwl(self):
        ret = Utils.check_gwl(self.files[0])
        self.assertIsNone(ret)
        
class Test_Map2Robot_main_single_barcode(unittest.TestCase):
    
    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.prefix = os.path.join(self.tmp_dir, 'output_')
        mapfile = os.path.join(data_dir, 'mapping_file_fecal_stability.txt')
        self.args = Map2Robot.parse_args(['--prefix', self.prefix, mapfile])
        self.files = Map2Robot.main(self.args)

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)        

    def test_main_gwl(self):
        ret = Utils.check_gwl(self.files[0])
        self.assertIsNone(ret)

        
if __name__ == '__main__':
    unittest.main()
