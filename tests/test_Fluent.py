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
from pyTecanFluent import Labware

# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


class Test_db(unittest.TestCase):
    def setUp(self):
        self.db = Fluent.db()

    def tearDown(self):
        pass

    def test_RackTypes(self):
        RackTypes = self.db.RackTypes()
        self.assertTrue(isinstance(RackTypes, list))

    def test_get_labware(self):
        RackType = self.db.RackTypes()[0]
        v = self.db.get_labware(RackType)
        self.assertTrue(isinstance(v, dict))

class Test_aspirate(unittest.TestCase):
    def setUp(self):
        self.asp = Fluent.aspirate()        

    def tearDown(self):
        pass

    def test_cmd(self):
        self.assertTrue(isinstance(self.asp.cmd, str))

class Test_dispense(unittest.TestCase):
    def setUp(self):
        self.disp = Fluent.dispense()

    def tearDown(self):
        pass

    def test_cmd(self):
        self.assertRaises(AssertionError, self.disp.cmd)

class Test_comment(unittest.TestCase):
    def setUp(self):
        self.comment = Fluent.comment('test')

    def tearDown(self):
        pass

    def test_cmd(self):
        ret = self.comment.cmd()
        self.assertTrue(isinstance(ret, str))
    
class Test_waste(unittest.TestCase):
    def setUp(self):
        self.waste = Fluent.waste()

    def tearDown(self):
        pass

    def test_cmd(self):
        ret = self.waste.cmd()
        self.assertTrue(isinstance(ret, str))        
            
class Test_gwl(unittest.TestCase):
    def setUp(self):
        self.gwl = Fluent.gwl()

    def tearDown(self):
        pass
    
    def test_add(self):
        asp = Fluent.aspirate()
        self.assertRaises(AssertionError, self.gwl.add, asp)
        
class Test_labware(unittest.TestCase):
    def setUp(self):
        self.labware = Labware.labware()

    def tearDown(self):
        pass

    def test_add(self):
        self.gwl = Fluent.gwl()
        ret = self.labware.add_gwl(self.gwl)
        self.assertTrue(ret is None)
    

if __name__ == '__main__':
    unittest.main()
