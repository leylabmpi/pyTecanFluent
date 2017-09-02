from __future__ import print_function

# import
import os
import sys
import numpy as np
import collections


# for adding labware, the following variables are needed
## labware_name = any name (must be unique)
## labware_type = one of the FluentControl labware definitions
## target_location = location type on worktable
## target_position = location position on worktable

# all labware definitions used 
LABWARE_DB = {
    'FCA, 50ul Filtered SBS' : {
        'target_location' : 'Nest7mm_Pos',
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

# boarder worktable sites that should NOT be used for asp/disp
BOARDER_SITES = [13,18,19,24,25,30,31,36]



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
        self._tip_cnt = collections.Counter()
        self._labware = collections.defaultdict(dict)
        self._tip_boxes = collections.defaultdict(dict)

    def add(self, x):
        """Adding labware/tip counts
        x: Fluent asp or disp object
        """
        # adding labware (labware_name : labware_type)
        if x.RackLabel in self._labware.keys():
            pass
        else:
            self._labware[x.RackLabel]['labware_type'] = x.RackType
            self._labware[x.RackLabel]['count'] = 1
            try: 
                self._labware[x.RackLabel]['target_location'] = LABWARE_DB[x.RackType]
            except KeyError:
                msg = 'WARNING: Labware type "{}" not recognized'
                print(msg.format(x.RackType), file=sys.stderr)
                self._labware[x.RackLabel]['target_location'] = None

        # adding tips
        if x.Volume < 45:
            tip_size = 50
        elif x.Volume < 180:
            tip_size = 200
        else:
            tip_size = 1000
        self._tip_cnt[self._tip_labware_type(tip_size)] += 1
        
    def _tip_labware_type(self, tip_size):
        return 'FCA, {}ul {}'.format(tip_size, self.tip_type)

    def _n_tip_boxes(self):
        """Counting number of tip boxes 
        """
        for labware_type,cnt in self._tip_cnt.items():
            box_cnt = int(round(cnt / 96.0 + 0.5, 0))    # assuming 96 tips per box
            for i in range(box_cnt):                
                labware_name = '{}[00{}]'.format(labware_type, i+1)                
                self._tip_boxes[labware_name]['labware_type'] = labware_type                
                self._tip_boxes[labware_name]['count'] = int(round(cnt / 96.0 + 0.5, 0))
                try: 
                    LABWARE_DB[labware_type]
                except KeyError:
                    msg = 'WARNING: Labware type "{}" not recognized'
                    print(msg.format(labware_type), file=sys.stderr)
                try: 
                    target_location = LABWARE_DB[labware_type]['target_location']
                except KeyError:
                    msg = 'WARNING: Labware type "{}" has no target location'
                    print(msg.format(labware_type), file=sys.stderr)
                    target_location = None
                self._tip_boxes[labware_name]['target_location'] = target_location            

    def labware_table(self):
        """making labware table with the following columns:
        labware_name
        labware_type
        labware_location
        labware_position
        """
        # counting tips & getting tip boxes
        self._n_tip_boxes()
        # adding tip boxes to labware table
        tbl = []
        site_cnt = 1
        for k in sorted(self._tip_boxes.keys()): 
            try: 
                labware_type = self._tip_boxes[k]['labware_type']
            except KeyError:
                labware_type = 'NA'
            try:
                target_location = self._tip_boxes[k]['target_location']
            except KeyError:
                target_location = 'NA'            
            tbl.append([k, labware_type, target_location, site_cnt])
            site_cnt += 1

        # adding other labware by target_position
        if site_cnt < 12:
            site_cnt = 12
        lw_tp = set([v['target_position'] for k,v in LABWARE_DB.items()])        
        for lw_cat in lw_cats:
            for k in sorted(self._labware.keys()):
                if self._labware[k]['category'] == lw_cat:
                    # site count
                    while 1:
                        site_cnt += 1                        
                        if site_cnt > max(BOARDER_SITES):
                            msg = 'Note enough worktable sites for all labware!'
                            raise ValueError
                        if site_cnt in BOURDER_SITES:
                            continue
                        else:
                            break                        
                    try: 
                        labware_type = self._tip_boxes[k]['labware_type']
                    except KeyError:
                        labware_type = 'NA'
                    try:
                        target_location = self._tip_boxes[k]['target_location']
                    except KeyError:
                        target_location = 'NA'            
                    tbl.append([k, labware_type, target_location, site_cnt])
                    site_cnt += 1
                                
        return tbl
        
    # @property
    # def Position(self):
    #     return self._Position

    # @Position.setter
    # def Position(self, value):
    #     self._Position = int(value)


# main
if __name__ == '__main__':
    pass
