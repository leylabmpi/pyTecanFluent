from __future__ import print_function

# import
import os
import sys
import numpy as np
import pandas as pd
import collections


# worktable-specific variables (move to config file?)
## all labware definitions used 
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

## total number of target_positions per target_location
TARGET_POSITIONS = {
    'Nest7mm_Pos': {
        'position_count' : 36, 
        'boarders' : (1,2,3,4,5,6,7,12,13,18,19,24,25,30,31,36)}, 
    'eppendorf' : {
        'position_count' : 96, 
        'boarders' : ()},
    'Trough_100ml_Wash_1' : {
        'position_count' : 8,
        'boarders' : ()}
}

class worktable_tracker():
    """Keeping track of available target positions for each target location
    """
    def __init__(self):
        self._target_positions = TARGET_POSITIONS
        # zero count for all avilable target_locations
        x = set([v['target_location'] for k,v in LABWARE_DB.items()])      
        x = {y:0 for y in x}
        self._tp_counter = collections.Counter(x)
        
    def add(self, target_location, n=1):
        """Adding item to worktable
        """
        try:
            self._tp_counter[target_location] += n
        except KeyError:
            msg = 'Target location "{}" not recognized'
            raise KeyError(msg.format(target_location))

    def next_available(self, target_location, add=True, boarder_check=True):
        """Getting next available target position for target location.
        add: add to target_position count
        boarder_check: screen out (skip) target_positions in the boarder
        """
        # stats on target_position for target_location
        i = self._current_count(target_location)
        max_pos = self._max_target_position(target_location)
        boarders = self._boarders(target_location)

        # finding next available position 
        while 1:
            i += 1    
            if i > max_pos:
                msg = 'Note enough target_positions for target_location: {}!'
                raise ValueError(msg.format(target_location))
            elif boarder_check == True and i in boarders:
                # skip boarder positions
                continue
            else:
                break
        # adding to count for target location 
        if add == True:
            self.add(target_location, 
                     i - self._current_count(target_location))
        # ret
        return i

    def _current_count(self, target_location):
        """Current count of target_positions filled for target location (1-indexing)
        """
        try:
            i = self._tp_counter[target_location]
        except KeyError:
            msg = 'Target location "{}" not recognized in target_position counter'
            raise KeyError(msg.format(target_location))
        return i
        
    def _max_target_position(self, target_location):
        """max target position for target location (1-indexing)
        """
        try:
            self._target_positions[target_location]
        except KeyError:
            msg = 'Target location "{}" not recognized in max target positions'
            raise KeyError(msg.format(target_location))
        try:
            max_pos = self._target_positions[target_location]['position_count']
        except KeyError:
            msg = 'Target location "{}" has no "position_count" key'
            raise KeyError(msg.format(target_location))
        return max_pos
        
    def _boarders(self, target_location):
        """Return board location (set) for target location
        """
        try:
            self._target_positions[target_location]
        except KeyError:
            msg = 'Target location "{}" not recognized in max target positions'
            raise KeyError(msg.format(target_location))
        try:
            boarders = self._target_positions[target_location]['boarders']
        except KeyError:
            msg = 'Target location "{}" has no "boarders" key'
            raise KeyError(msg.format(target_location))
        return boarders

    def all_target_locations(self):
        """Return all target_positions
        """
        return self._tp_counter.keys()


    
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

    def add(self, x, add_tip=True):
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
                self._labware[x.RackLabel]['target_location'] = \
                    LABWARE_DB[x.RackType]['target_location']
            except KeyError:
                msg = 'WARNING: Labware type "{}" not recognized'
                print(msg.format(x.RackType), file=sys.stderr)
                self._labware[x.RackLabel]['target_location'] = None

        if add_tip == True:
            # adding tips
            if x.Volume < 45:
                tip_size = 50
            elif x.Volume < 180:
                tip_size = 200
            else:
                tip_size = 1000
            self._tip_cnt[self._tip_labware_type(tip_size)] += 1
        
    def _tip_labware_type(self, tip_size):
        """Getting tip box labware type
        tip_size: size of tip (in ul)
        """
        return 'FCA, {}ul {}'.format(tip_size, self.tip_type)

    def _n_tip_boxes(self):
        """Counting number of tip boxes and adding tips boxes to labware
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
        # counting tips & adding tip boxes to labware
        self._n_tip_boxes()
        # adding labware by target location 
        wt_tracker = worktable_tracker()

        # addig tip boxes to labware table
        tbl = []
        for labware_name in sorted(self._tip_boxes.keys()):             
            # getting info for box with labware name <k>
            try: 
                labware_type = self._tip_boxes[labware_name]['labware_type']
            except KeyError:
                labware_type = 'NA'
            try:
                target_location = self._tip_boxes[labware_name]['target_location']
            except KeyError:
                target_location = 'NA'           
            # getting target_position
            target_position = wt_tracker.next_available(target_location, 
                                                        boarder_check=False)
            tbl.append([labware_name, labware_type, 
                        target_location, target_position])
        
        # adding other labware by target_position
        for target_location in wt_tracker.all_target_locations():
            # all labware to add (besides tips)  
            for labware_name in self._labware.keys():
                if self._labware[labware_name]['target_location'] == target_location:
                    # labware type
                    try: 
                        labware_type = self._labware[labware_name]['labware_type']
                    except KeyError:
                        labware_type = 'NA'
                    # target_position
                    target_position = wt_tracker.next_available(target_location)
                    # adding to labware table
                    tbl.append([labware_name, labware_type, 
                                target_location, target_position])

        # convert to pandas data.frame
        labels = ['labware_name', 'labware_type', 'target_location', 'target_position']
        df = pd.DataFrame.from_records(tbl, columns=labels)
                                
        return df
        

# main
if __name__ == '__main__':
    pass
