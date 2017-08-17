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
from pyTecanFluent import Dilute
from pyTecanFluent import Utils


# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
class Test_Dilute_import(unittest.TestCase):

    def setUp(self):
        concfile = os.path.join(data_dir, 'conc_file1.txt')
        args = Dilute.parse_args([concfile])
        Dilute.check_args(args)
        self.df_conc = Dilute.conc2df(args.concfile)

    def tearDown(self):
        self.df_conc = None
        pass

    # import
    def test_load_conc_txt(self):
        self.assertTrue(isinstance(self.df_conc, pd.DataFrame))


class Test_Dilute_rows(unittest.TestCase):

    def setUp(self):
        concfile = os.path.join(data_dir, 'conc_file1.txt')
        self.args = Dilute.parse_args(['--rows', '1,2,4-5', concfile])
        Dilute.check_args(self.args)
        self.df_conc = Dilute.conc2df(self.args.concfile,
                                     row_select=self.args.rows)

    def tearDown(self):
        self.df_conc = None
        pass

    # import
    def test_load_conc_txt(self):
        self.assertTrue(isinstance(self.df_conc, pd.DataFrame))
        self.assertEqual(self.df_conc.shape[0], 4)
        self.assertListEqual(list(self.df_conc.index), [0,1,3,4])

class Test_Dilute_check_df_conc(unittest.TestCase):

    def setUp(self):
        concfile = os.path.join(data_dir, 'conc_file1.txt')
        self.args = Dilute.parse_args([concfile])
        Dilute.check_args(self.args)
        self.df_conc = Dilute.conc2df(self.args.concfile)

    def tearDown(self):
        self.df_conc = None
        pass

    # import
    def test_load_conc_txt(self):
        self.assertTrue(isinstance(self.df_conc, pd.DataFrame))
        Dilute.check_df_conc(self.df_conc, self.args)

class Test_Dilute_dil_vols(unittest.TestCase):

    def setUp(self):
        concfile = os.path.join(data_dir, 'conc_file1.txt')
        self.args = Dilute.parse_args([concfile])
        Dilute.check_args(self.args)
        self.df_conc = Dilute.conc2df(self.args.concfile)
        Dilute.check_df_conc(self.df_conc, self.args)
        self.df_conc = Dilute.dilution_volumes(self.df_conc,
                                               dilute_conc=self.args.dilution,
                                               min_vol=self.args.minvolume,
                                               max_vol=self.args.maxvolume,
                                               min_total=self.args.mintotal,
                                               dest_type=self.args.desttype)

    def tearDown(self):
        self.df_conc = None
        pass

    # import
    def test_load_conc_txt(self):
        self.assertTrue(isinstance(self.df_conc, pd.DataFrame))


class Test_Dilute_addDest(unittest.TestCase):

    def setUp(self):
        concfile = os.path.join(data_dir, 'conc_file1.txt')
        self.args = Dilute.parse_args([concfile])
        Dilute.check_args(self.args)
        self.df_conc = Dilute.conc2df(self.args.concfile)

    def tearDown(self):
        self.df_conc = None
        self.args = None

    # adding destination
    def test_load_conc_txt(self):
        self.df_conc = Dilute.add_dest(self.df_conc, self.args.dest, self.args.deststart)
        self.assertTrue(isinstance(self.df_conc, pd.DataFrame))

    # destination start offset
    def test_load_conc_deststart(self):
        self.df_conc = Dilute.add_dest(self.df_conc,
                                       self.args.dest,
                                       dest_start=49)

        self.assertTrue(isinstance(self.df_conc, pd.DataFrame))
        loc_start = self.df_conc.loc[0,'dest_location']
        self.assertEqual(loc_start, 49)


class Test_Dilute_main1(unittest.TestCase):

    def setUp(self):
        concfile = os.path.join(data_dir, 'conc_file1.txt')
        self.args = Dilute.parse_args(['--prefix', '/tmp/DIL1', concfile])
        self.files = Dilute.main(self.args)

    def tearDown(self):
        pass

    def test_main_gwl(self):
        ret = Utils.check_gwl(self.files[0])
        self.assertIsNone(ret)


class Test_Dilute_main2(unittest.TestCase):

    def setUp(self):
        concfile = os.path.join(data_dir, 'conc_file2.txt')
        self.args = Dilute.parse_args(['--prefix', '/tmp/DIL2', concfile])
        self.files = Dilute.main(self.args)

    def tearDown(self):
        pass

    def test_main_gwl(self):
        ret = Utils.check_gwl(self.files[0])
        self.assertIsNone(ret)
