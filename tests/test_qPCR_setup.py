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
from pyTecanFluent import QPCR
from pyTecanFluent import Utils


# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')



# tests
class Test_import(unittest.TestCase):

    def setUp(self):
        infile = os.path.join(data_dir, 'qPCR_setup/qPCR_Zach_plate1.xlsx')
        args = QPCR.parse_args([infile])
        QPCR.check_args(args)
        self.df_setup = QPCR.load_setup(args.setup)
        QPCR.check_df_setup(self.df_setup)

    def tearDown(self):
        self.df_setup = None
        pass

    # import
    def test_load_excel(self):
        self.assertTrue(isinstance(self.df_setup, pd.DataFrame))
        #print(self.df_setup)


class Test_add_dest(unittest.TestCase):

    def setUp(self):
        infile = os.path.join(data_dir, 'qPCR_setup/qPCR_Zach_plate1.xlsx')
        args = QPCR.parse_args([infile])
        QPCR.check_args(args)
        self.df_setup = QPCR.load_setup(args.setup)
        QPCR.check_df_setup(self.df_setup)
        QPCR.add_dest(self.df_setup, args.dest, dest_type=args.desttype)
        
    def tearDown(self):
        self.df_setup = None
        pass

    # import
    def test_add_dest(self):
        self.assertTrue(isinstance(self.df_setup, pd.DataFrame))
        #print(self.df_setup)

class Test_qPCR_main1(unittest.TestCase):

    def setUp(self):
        infile = os.path.join(data_dir, 'qPCR_setup/qPCR_Zach_plate1.xlsx')
        args = QPCR.parse_args(['--prefix', '/tmp/qPCR', infile])
        self.files = QPCR.main(args)

    def tearDown(self):
        pass

    def test_main_gwl(self):
        ret = Utils.check_gwl(self.files[0])
        self.assertIsNone(ret)
