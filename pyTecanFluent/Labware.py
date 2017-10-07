from __future__ import print_function

# import
import os
import sys
import string
import itertools
import numpy as np
import pandas as pd
import collections


# worktable-specific variables (move to config file?)
## all labware definitions used 
LABWARE_DB = {
    'FCA, 10ul SBS High' : {
        'target_location' : 'Nest7mm_Pos',
        'category' : 'tip',
        'wells' : 96,
        'max_volume' : None},
    'FCA, 50ul SBS High' : {
        'target_location' : 'Nest7mm_Pos',
        'category' : 'tip',
        'wells' : 96,
        'max_volume' : None},
    'FCA, 200ul SBS High' : {
        'target_location' : 'Nest7mm_Pos',
        'category' : 'tip',
        'wells' : 96,
        'max_volume' : None},
    'FCA, 1000ul SBS High' : {
        'target_location' : 'Nest7mm_Pos',
        'category' : 'tip',
        'wells' : 96,
        'max_volume' : None},
    'FCA, 10ul Filtered SBS High' : {
        'target_location' : 'Nest7mm_Pos',
        'category' : 'tip',
        'wells' : 96,
        'max_volume' : None},
    'FCA, 50ul Filtered SBS High' : {
        'target_location' : 'Nest7mm_Pos',
        'category' : 'tip',
        'wells' : 96,
        'max_volume' : None},
    'FCA, 200ul Filtered SBS High' : {
        'target_location' : 'Nest7mm_Pos',
        'category' : 'tip',
        'wells' : 96,
        'max_volume' : None},
    'FCA, 1000ul Filtered SBS High' : {
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
        'boarders' : (1,2,3,4,5,6,7,8,9,10,11,12,
                      13,18,19,24,25,30,31,36)}, 
    'eppendorf' : {
        'position_count' : 96, 
        'boarders' : ()},
    'Trough_100ml_Wash_1' : {
        'position_count' : 8,
        'boarders' : ()}
}

def total_positions(labware_type):
    """Getting total number of positions for the labware type
    """
    try:
        LABWARE_DB[labware_type]
    except KeyError:
        msg = 'Labware type "{}" not recognised'
        raise KeyError(msg.format(labware_type))
    try:
        positions = LABWARE_DB[labware_type]['wells']
    except KeyError:
        msg = 'No "wells" key for labware type "{}"'
        raise KeyError(msg.format(labware_type))
    return positions

def position2well(position, wells=96, just_row=False, just_col=False):
    """Convert position to well
    Note: assuming column-wise ordering
    """
    # making plate index
    wells = int(wells)
    if wells == 96:
        nrows = 8
        ncols = 12
    elif wells == 384:
        nrows = 16
        ncols = 24
    else:
        raise ValueError('Number of wells ({}) not recognized'.format(wells))
    
    rows = list(string.ascii_uppercase[:nrows])
    cols = [x + 1 for x in range(ncols)]
    # position : [row, col]
    pos_idx = {i+1:x for i,x in enumerate(itertools.product(cols, rows))}
    # getting row-col
    try:
        col,row = pos_idx[position]
    except KeyError:
        msg = 'Cannot find well for position: {}'
        raise KeyError(msg.format(position))
    if just_row == True:
        return row
    elif just_col == True:
        return col
    else:
        return [row,col]


def well2position(well, wells=96, labware_type=None):
    """well ID to column-wise position
    """
    # setting wells based on labware type (if provided)
    if labware_type is not None:
        wells = total_positions(labware_type)
    # skip if already integer (not well)
    try:
        return int(well)
    except ValueError:
        pass
    # well --> position index
    wells = int(wells)
    if wells == 96:
        nrows = 8
        ncols = 12
    elif wells == 384:
        nrows = 16
        ncols = 24
    else:
        raise ValueError('Number of wells ({}) not recognized'.format(wells))    
    rows = list(string.ascii_uppercase[:nrows])
    cols = [x + 1 for x in range(ncols)]
    well_idx = {'{0}{1:0>2}'.format(x[1], x[0]):i+1 for i,x in enumerate(itertools.product(cols, rows))}
    # getting position
    try:
        position = well_idx[well]
    except KeyError:
        msg = 'Cannot find well "{}"'
        raise KeyError(msg.format(well))
    # return 
    return position


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
        target_location : int, target location
        n : int, how much to add
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
            self.add(target_location, i - self._current_count(target_location))
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
    def __init__(self, tip_types={}):
        self.tip_types = tip_types
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
        # adding tip
        if add_tip == True:
            # adding tips
            if x.Volume < 45:
                tip_size = 50
            elif x.Volume < 180:
                tip_size = 200
            else:
                tip_size = 1000
            
            self._tip_cnt[self._tip_labware_type(tip_size)] += 1

    def next_target_position(self, labware_type, boarder_check=False):
        """Getting the next target position for the labware type
        """
        # getting sum for labware type
        i = len([k for k,v in self._labware.items() if v['labware_type'] == labware_type])
        try:
            target_location = LABWARE_DB[labware_type]['target_location']
        except KeyError:
            msg = 'Labware type "{}" not recognized'
            raise KeyError(msg.format(labware_type))
        
        wt_tracker = worktable_tracker()        
        max_pos = wt_tracker._max_target_position(target_location)
        # finding next available position
        if boarder_check == True:
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
        else:
            i += 1
            if i > max_pos:
                msg = 'Note enough target_positions for target_location: {}!'
                raise ValueError(msg.format(target_location))
        return i
                        
    def _tip_labware_type(self, tip_size):
        """Getting tip box labware type
        tip_size: size of tip (in ul)
        """
        try:
            return self.tip_types[tip_size]
        except KeyError:
            msg = 'ERROR: no tip type for tip size: {}'
            raise KeyError(msg.format(tip_size))                       

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
