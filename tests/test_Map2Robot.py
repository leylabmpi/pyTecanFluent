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
from pyTecanFluent import Map2Robot
from pyTecanFluent import Utils


# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
class Test_Map2Robot_parser(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    # args
    def test_parser_dest(self):
        # 96-well
        args = Map2Robot.parse_args(['--desttype', '96 Well Eppendorf TwinTec PCR',
                                     '--deststart', '97',
                                     'mapfile'])
        with self.assertRaises(ValueError):
            Map2Robot.check_args(args)
        # 384-well
        args = Map2Robot.parse_args(['--desttype', '384 Well Biorad PCR',
                                     '--deststart', '389',
                                     'mapfile'])
        with self.assertRaises(ValueError):
            Map2Robot.check_args(args)

    def test_parser_tube(self):
        # mastermix <1
        args = Map2Robot.parse_args(['--mmtube', '0', 'mapfile'])
        with self.assertRaises(ValueError):
            Map2Robot.check_args(args)
        # mastermix >24
        args = Map2Robot.parse_args(['--mmtube', '25', 'mapfile'])
        with self.assertRaises(ValueError):
            Map2Robot.check_args(args)

class Test_Map2Robot_import(unittest.TestCase):

    def setUp(self):
        mapfile = os.path.join(data_dir, 'mapping_file_fecal_stability.txt')
        args = Map2Robot.parse_args([mapfile])
        Map2Robot.check_args(args)
        self.df_map = Map2Robot.map2df(args.mapfile)

    def tearDown(self):
        self.df_map = None
        pass

    # import
    def test_load_map_txt(self):
        self.assertTrue(isinstance(self.df_map, pd.DataFrame))

    def test_load_map_xls(self):
        self.assertTrue(isinstance(self.df_map, pd.DataFrame))

class Test_Map2Robot_rows(unittest.TestCase):

    def setUp(self):
        mapfile = os.path.join(data_dir, 'mapping_file_fecal_stability.txt')
        self.args = Map2Robot.parse_args(['--rows', '1,2,4-5', mapfile])
        Map2Robot.check_args(self.args)
        self.df_map = Map2Robot.map2df(self.args.mapfile,
                                       row_select=self.args.rows)

    def tearDown(self):
        self.df_map = None
        pass

    # import
    def test_load_map_txt(self):
        self.assertTrue(isinstance(self.df_map, pd.DataFrame))
        self.assertEqual(self.df_map.shape[0], 4)
        self.assertListEqual(list(self.df_map.index), [0,1,3,4])

class Test_Map2Robot_addDest(unittest.TestCase):

    def setUp(self):
        mapfile = os.path.join(data_dir, 'mapping_file_fecal_stability.txt')
        self.args = Map2Robot.parse_args([mapfile])
        Map2Robot.check_args(self.args)
        self.df_map = Map2Robot.map2df(self.args.mapfile)

    def tearDown(self):
        self.df_map = None
        self.args = None

    # adding destination
    def test_load_map_txt(self):
        self.df_map = Map2Robot.add_dest(self.df_map,
                                         self.args.destname)
        self.assertTrue(isinstance(self.df_map, pd.DataFrame))

    # destination start offset
    def test_load_map_deststart(self):
        self.df_map = Map2Robot.add_dest(self.df_map,
                                         self.args.destname,
                                         dest_start=49)

        self.assertTrue(isinstance(self.df_map, pd.DataFrame))
        loc_start = self.df_map.loc[0,'TECAN_dest_target_position']
        self.assertEqual(loc_start, 49)

class Test_Map2Robot_destStart384(unittest.TestCase):

    def setUp(self):
        mapfile = os.path.join(data_dir, 'mapping_file_fecal_stability.txt')
        self.args = Map2Robot.parse_args(['--desttype',
                                          '384 Well Biorad PCR',
                                          mapfile])
        Map2Robot.check_args(self.args)
        self.df_map = Map2Robot.map2df(self.args.mapfile)

    def tearDown(self):
        self.df_map = None
        self.args = None

    def test_load_map_deststart384(self):
        self.df_map = Map2Robot.add_dest(self.df_map,
                                         dest_labware=self.args.destname,
                                         dest_type=self.args.desttype,
                                         dest_start=200)
        self.assertTrue(isinstance(self.df_map, pd.DataFrame))
        loc_start = self.df_map.loc[0,'TECAN_dest_target_position']
        self.assertEqual(loc_start, 200.0)

    def test_load_map_deststart384_2(self):
        self.df_map = Map2Robot.add_dest(self.df_map,
                                         dest_labware=self.args.destname,
                                         dest_type=self.args.desttype,
                                         dest_start=370)
        self.assertTrue(isinstance(self.df_map, pd.DataFrame))

        loc_start = self.df_map.loc[0,'TECAN_dest_target_position']
        self.assertEqual(loc_start, 370.0)
        
        #i = self.df_map.shape[0]-1
        #loc_end = self.df_map.ix[i,'TECAN_dest_target_position']
        #self.assertEqual(loc_end, 384.0)

    # destination location reorder for 384-well plates
    def test_load_map_destreorder(self):
        self.df_map = Map2Robot.add_dest(self.df_map,
                                         dest_labware=self.args.destname)
        self.df_map = Map2Robot.reorder_384well(self.df_map,
                                                'TECAN_dest_target_position')

        self.assertTrue(isinstance(self.df_map, pd.DataFrame))
        loc_start = self.df_map.loc[0,'TECAN_dest_target_position']
        self.assertEqual(loc_start, 1)
        nrow = self.df_map.shape[0]
        loc_end = self.df_map.loc[nrow-1,'TECAN_dest_target_position']
        self.assertEqual(loc_end, 96)

class Test_Map2Robot_pipMasterMix(unittest.TestCase):

    def setUp(self):
        mapfile = os.path.join(data_dir, 'mapping_file_fecal_stability.txt')
        self.args = Map2Robot.parse_args([mapfile])
        Map2Robot.check_args(self.args)
        self.df_map = Map2Robot.map2df(self.args.mapfile)
        self.df_map = Map2Robot.add_dest(self.df_map,
                                         dest_labware=self.args.destname,
                                         dest_type=self.args.desttype)
        if self.args.desttype == '384':
            self.df_map = reorder_384well(self.df_map, 'TECAN_dest_target_position')
        self.gwl = Fluent.gwl()
        self.gwl_file = '/tmp/MM.gwl'

    def tearDown(self):
        self.df_map = None
        self.args = None

    def test_pip_master_mix(self):
        Map2Robot.pip_mastermix(self.df_map, self.gwl)
        self.gwl.write(self.gwl_file)
        ret = Utils.check_gwl(self.gwl_file)
        self.assertIsNone(ret)

class Test_Map2Robot_main_basic(unittest.TestCase):

    def setUp(self):
        mapfile = os.path.join(data_dir, 'basic_96well.txt')
        self.args = Map2Robot.parse_args(['--prefix', '/tmp/MAP', mapfile])
        self.files = Map2Robot.main(self.args)

    def tearDown(self):
        pass

    def test_main_gwl(self):
        ret = Utils.check_gwl(self.files[0])
        self.assertIsNone(ret)
        
class Test_Map2Robot_main_single_barcode(unittest.TestCase):

    def setUp(self):
        mapfile = os.path.join(data_dir, 'mapping_file_fecal_stability.txt')
        self.args = Map2Robot.parse_args(['--prefix', '/tmp/MAP', mapfile])
        self.files = Map2Robot.main(self.args)

    def tearDown(self):
        pass

    def test_main_gwl(self):
        ret = Utils.check_gwl(self.files[0])
        self.assertIsNone(ret)

        
if __name__ == '__main__':
    unittest.main()
