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
class Test_Pool_import(unittest.TestCase):

    def setUp(self):
        self.pcr_file = os.path.join(data_dir, 'PCR-run1.xlsx')
        self.map_file = os.path.join(data_dir, 'samples_S5-N7.xlsx')
        self.args = Pool.parse_args([self.pcr_file])
        Pool.check_args(self.args)

    def tearDown(self):
        self.df_conc = None

    # import
    def test_load_sample_files(self):
        df_samps = []
        for f in self.args.samplefiles:
            df_samp = Pool.sample2df(f,
                                     sample_col=self.args.sample_col,
                                     include_col=self.args.include_col,
                                     labware_name_col=self.args.sample_labware_name,
                                     labware_type_col=self.args.sample_labware_type,
                                     position_col=self.args.position_col,
                                     file_format=self.args.sample_format,
                                     row_select=self.args.sample_rows, 
                                     header=self.args.sample_header)
            df_samps.append(df_samp)
        df_samp = pd.concat(df_samps)
        self.assertTrue(isinstance(df_samp, pd.DataFrame))

    def test_load_mapping_file(self):
        df_map = Pool.map2df(self.map_file,                    
                             file_format=self.args.map_format,
                             header=self.args.map_header)
        self.assertTrue(isinstance(df_map, pd.DataFrame))


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
