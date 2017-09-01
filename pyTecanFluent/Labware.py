from __future__ import print_function

# import
import os
import sys
import numpy as np
import collections
#from collections import Counter


# for adding labware, the following variables are needed
## labware_name = any name (must be unique)
## labware_type = one of the FluentControl labware definitions
## target_location = location type on worktable
## target_position = location position on worktable

LABWARE_DB = {
    'FCA, 50ul Filtered SBS' : {
        'labware_type' : 'Nest7mm_Pos',
        'category' : 'tip',
        'wells' : 96,
        'max_volume' : None},
    'FCA, 200ul Filtered SBS' : {
        'target_location' : 'Nest7mm_Pos',
        'category' : 'tip',
        'wells' : 96,
        'max_volume' : None},
    'FCA, 1000ul Filtered SBS' : {
        'target_location' : 'Nest7mm_Pos',
        'category' : 'tip',
        'wells' : 96,
        'max_volume' : None},
    '100ml_1' : {
        'target_location' : 'Trough_100ml_Wash_1',
        'category' : 'trough',
        'wells' : 1,
        'max_volume' : 100000},
    '96 Well Eppendorf TwinTec PCR' : {
        'target_location' : 'Nest7mm_Pos',
        'category' : 'plate',
        'wells' : 96,
        'max_volume' : 200},
    '384 Well Biorad PCR' : {
        'target_location' : 'Nest7mm_Pos',
        'category' : 'plate',
        'wells' : 384,
        'max_volume' : 40},
    'OptiPlate_96F' : {
        'target_location' : 'Nest7mm_Pos',
        'category' : 'plate',
        'wells' : 96,
        'max_volume' : 300},
    '1.5ml Eppendorf' : {
        'target_location' : 'eppendorf',
        'category' : 'tube',
        'wells' : 1,
        'max_volume' : 145000},
    '2.0ml Eppendorf' : {
        'target_location' : 'eppendorf',
        'category' : 'tube',
        'wells' : 1,
        'max_volume' : 195000}
}



class labware_tracker():
    """Class for tracking which labware is needed for the worklist commands.
    The number of tips are counted.
    The labware is counted.

    Needed variables for adding labware to worktable:
    *) labware_name = any name (must be unique) 
    *) labware_type = one of the FluentControl labware definitions
    *) target_location = location type on worktable
    *) target_position = location position on worktable
    """
    def __init__(self):
        self.tip_type = 'Filtered SBS'
        self.tip_cnt = collections.Counter()
        self.labware = collections.defaultdict(dict)

    def add(self, x):
        """Adding labware/tip counts
        x: Fluent asp or disp object
        """
        # adding labware
        ## needed: 
        if x.RackLabel in self.labware.keys():
            pass
        else:
            self.labware[x.RackLabel]['RackType'] = x.RackType
        # adding tips
        if x.Volume < 45:
            tip_size = 50
        elif x.Volume < 180:
            tip_size = 200
        else:
            tip_size = 1000
        self.tip_cnt[self._tip_labware_type(tip_size)] += 1
        
    def _tip_labware_type(self, tip_size):
        return '{}ul {}'.format(tip_size, self.tip_type)

    # TODO: make labware output (with header)

# main
if __name__ == '__main__':
    pass
