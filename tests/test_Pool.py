#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import
## batteries
import os
import sys
import unittest
## 3rd party
import pandas as pd
## package
from pyTecanFluent import Pool


# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
class Test_96well_NoMap(unittest.TestCase):

    def setUp(self):
        self.pcr_file = os.path.join(data_dir, 'PCR-run1.xlsx')
        self.args = Pool.parse_args(['--prefix', '/tmp/POOL1',
                                     self.pcr_file])

    def tearDown(self):
        self.df_conc = None

    def test_main(self):
        Pool.main(self.args)

class Test_2_96well_Map(unittest.TestCase):

    def setUp(self):
        self.pcr_file1 = os.path.join(data_dir, 'PCR-run1.xlsx')
        self.pcr_file2 = os.path.join(data_dir, 'PCR-run2.xlsx')
        self.map_file = os.path.join(data_dir, 'samples_S5-N7.xlsx')
        self.args = Pool.parse_args(['--prefix', '/tmp/POOL2',
                                     '--mapfile', self.map_file,
                                     self.pcr_file1, self.pcr_file2])

    def tearDown(self):
        self.df_conc = None

    def test_main(self):
        Pool.main(self.args)


if __name__ == '__main__':
    unittest.main()
