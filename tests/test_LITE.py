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
from pyTecanFluent import LITE
from pyTecanFluent import Utils


# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests        
class Test_LITE_main(unittest.TestCase):

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.prefix = os.path.join(self.tmp_dir, 'output_')
        concfile = os.path.join(data_dir, 'basic_96well.txt')
        self.args = LITE.parse_args(['--prefix', self.prefix, concfile])
        self.files = LITE.main(self.args)

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)               

    def test_main_gwl(self):
        ret = Utils.check_gwl(self.files[0])
        self.assertIsNone(ret)

class Test_LITE_onePlate(unittest.TestCase):

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.prefix = os.path.join(self.tmp_dir, 'output_')
        concfile = os.path.join(data_dir, 'LITE_1plate.xlsx')
        self.args = LITE.parse_args(['--prefix', self.prefix,
                                     '--pcr-mm-volume', '7', 
                                     '--pcr-mm-labware-type', '10ml Falcon',
                                     concfile])
        self.files = LITE.main(self.args)

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)               

    def test_main_gwl(self):
        ret = Utils.check_gwl(self.files[0])
        self.assertIsNone(ret)

if __name__ == '__main__':
    unittest.main()
