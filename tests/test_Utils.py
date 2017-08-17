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
from pyTecanFluent import Utils


# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
class Test_Utils(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_make_range(self):
        x = Utils.make_range('all')
        self.assertIsNone(x)

        x = Utils.make_range('0')
        self.assertListEqual(x, [0])

        x = Utils.make_range('1,2,5')
        self.assertListEqual(x, [1,2,5])

        x = Utils.make_range('1,2,5-6')
        self.assertListEqual(x, [1,2,5,6])

        x = Utils.make_range('1-3,5-6')
        self.assertListEqual(x, [1,2,3,5,6])


    def test_make_range_zeroindex(self):
        x = Utils.make_range('all', set_zero_index=True)
        self.assertIsNone(x)

        with self.assertRaises(ValueError):
            Utils.make_range('0', set_zero_index=True)

        x = Utils.make_range('1,2,5', set_zero_index=True)
        self.assertListEqual(x, [0,1,4])

        x = Utils.make_range('1,2,5-6', set_zero_index=True)
        self.assertListEqual(x, [0,1,4,5])

    def test_check_gwl(self):
        gwl_file = os.path.join(data_dir, 'multi_dispense.gwl')
        ret = Utils.check_gwl(gwl_file)
        self.assertIsNone(ret)
