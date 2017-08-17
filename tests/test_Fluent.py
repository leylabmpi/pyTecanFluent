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
from pyTecanFluent import Fluent


# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
class Test_Fluent(unittest.TestCase):

    def setUp(self):
        self.rd = Fluent.reagent_distribution()
        self.asp = Fluent.aspirate()
        self.disp = Fluent.dispense()

    def tearDown(self):
        pass

    def test_init_rd(self):
        cmd = self.rd.cmd()
        self.assertTrue(isinstance(cmd, str))

    def test_init_asp(self):
        cmd = self.asp.cmd()
        self.assertTrue(isinstance(cmd, str))

    def test_init_rd(self):
        cmd = self.disp.cmd()
        self.assertTrue(isinstance(cmd, str))
