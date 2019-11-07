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
from pyTecanFluent import Dilute
from pyTecanFluent import Utils


# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
class Test_Dilute_main1(unittest.TestCase):

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.prefix = os.path.join(self.tmp_dir, 'output_')
        concfile = os.path.join(data_dir, 'conc_file1.txt')
        self.args = Dilute.parse_args(['--prefix', self.prefix, concfile])
        self.files = Dilute.main(self.args)

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)        

    def test_main_gwl(self):
        ret = Utils.check_gwl(self.files[0])
        self.assertIsNone(ret)
        
class Test_Dilute_main2(unittest.TestCase):

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.prefix = os.path.join(self.tmp_dir, 'output_')
        concfile = os.path.join(data_dir, 'conc_file2.txt')
        self.args = Dilute.parse_args(['--prefix', '/tmp/DIL2', concfile])
        self.files = Dilute.main(self.args)

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)        

    def test_main_gwl(self):
        ret = Utils.check_gwl(self.files[0])
        self.assertIsNone(ret)


if __name__ == '__main__':
    unittest.main()
