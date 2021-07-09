from __future__ import print_function

# import
## batteries
import os
import sys
import json
import string
import itertools
import collections
## 3rd party
import numpy as np
import pandas as pd
## package
from pyTecanFluent import Fluent


class utils(object):
    """Utility functions for labware
    """

    def __init__(self):
        # target position
        d = os.path.join(os.path.split(__file__)[0],  'database')
        f = os.path.join(d, 'labware.json')
        with open(f) as inF:
            self.labware = json.load(inF)

    def get_wells(self, RackType):
        """Getting wells of RackType
        """
        if RackType is None:
            return None
        try:
            wells = self.labware[RackType]['wells']
        except KeyError:
            return None
        return wells
            
    def position2well(self, position, wells=96, just_row=False, just_col=False):
        """Convert position to well
        Note: assuming column-wise ordering
        Return: list with row & column IDs => [row,column]
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
            msg = 'Number of wells ({}) not recognized'
            raise ValueError(msg.format(wells))
    
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

    def well2position(self, well, wells=96, RackType=None):
        """well ID to column-wise position (opposite of position2well)
        """
        # skip if already integer (not wellID)
        try:
            return int(well)
        except ValueError:
            pass
        # if RackType provided, selecting wells from database
        if RackType is not None:
            try:
                wells = self.labware[RackType]['wells']
            except KeyError:
                msg = 'No wells found for RackType: "{}"'
                raise KeyError(msg.format(RackType))
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
            try:
                # adding zero-padding
                position = well_idx[well[0] + '0' + well[1:]]
            except KeyError:
                msg = 'Cannot find well "{}"'
                raise KeyError(msg.format(well))
        # return 
        return position
    
class labware(object):
    """Class for summarizing labware in a gwl object.
    Note: tip boxes are considered separate from other labware
    """
    def __init__(self):
        self.tip_count = {}
        self.tip_boxes = {}
        self.labware = {} 
        self.labware_order = {}
        # target position
        d = os.path.join(os.path.split(__file__)[0],  'database')
        f = os.path.join(d, 'target_position.json')
        with open(f) as inF:
            self.target_position = json.load(inF)
                
    def add_gwl(self, gwl):
        """Adding labware from gwl object to labware object.
        Note: this can be used to sum up labware from multiple gwl objects.
        """
        # counting tips
        self._count_tips(gwl.commands)
        # adding labware 
        for cmd in gwl.commands:
            self._add_labware(cmd, gwl)
        # summing up tip boxes        
        self._add_tip_boxes(gwl)

    def table(self):
        """Creating pandas dataframe of labware
        columns: labware_name, labware_type,target_location,target_position
        """
        # init target position counters
        loc_tracker = {}
        
        # init liist of dicts (will be coverted to dataframe)
        cols = ['labware_name', 'labware_type',
                'target_location', 'target_position']
        
        # adding tip boxes; sorting largest to smallest tip size; [001], [002], ...
        df_tips = []
        func = lambda x: (x[1][1]['max_volume'], x[0][0])
        for RackLabel,v in sorted(self.tip_boxes.items(), key=func, reverse=True):
            RackType_num,v = v
            # RackType
            try:
                RackType = v['RackType']
            except KeyError:
                msg = 'No RackType for labware: "{}"'
                raise KeyError(msg.format(RackLabel))
            # location & position
            loc,pos = self._next(v, loc_tracker, keep_empty=True)
            if loc is None:
                msg = 'No possible target location for labware: "{}"'
                raise ValueError(msg.format(RackLabel))
            # creating table entry
            df_tips.append({'labware_name' : RackLabel,
                            'labware_type' : RackType,
                            'target_location' : loc,
                            'target_position' : pos,
                            'target_location_prompt' : loc,
                            'target_position_prompt' : pos})  # prompt = what is prompted for user
        # adding other labware; sorting by labware order in the gwl
        df_labware = []
        adapter_cnt = {'96 well' : 0, '384 well' : 0, '96 well magnet' : 0}
        for RackLabel in sorted(self.labware_order, key=self.labware_order.get):
            v = self.labware[RackLabel]
            # RackType
            try:
                RackType = v['RackType']
            except KeyError:
                msg = 'No RackType for labware: "{}"'
                raise KeyError(msg.format(RackLabel))
            # location & position
            loc,pos = self._next(v, loc_tracker, keep_empty=True)
            loc_prompt = loc
            pos_prompt = pos
            if loc is None:
                msg = 'No possible target location for labware: "{}"'
                raise ValueError(msg.format(RackLabel))
            ## if Racktype includes adapter, adding adapter + plate
            if RackType.startswith('PCR Adapter 96 Well and '):
                adapter_cnt['96 well'] += 1
                # adding adapter to labware table
                df_labware.append({'labware_name' : 'PCR Adapter for {}'.format(RackLabel),
                                   'labware_type' : 'PCR Adapter 96 Well',
                                   'target_location' : loc,
                                   'target_position' : pos,
                                   'target_location_prompt' : loc,
                                   'target_position_prompt' : pos})
                # editing info for plate on adapter
                RackType = RackType.replace('PCR Adapter 96 Well and ', '')
                loc = 'PCR96WellAdapter_CoverSite'
                pos = adapter_cnt['96 well']
            elif RackType.startswith('PCR Adapter 384 Well and '):
                adapter_cnt['384 well'] += 1
                # adding adapter to labware table
                df_labware.append({'labware_name' : 'PCR Adapter for {}'.format(RackLabel),
                                   'labware_type' : 'PCR Adapter 384 Well',
                                   'target_location' : loc,
                                   'target_position' : pos,
                                   'target_location_prompt' : loc,
                                   'target_position_prompt' : pos})
                # editing info for plate on adapter
                RackType = RackType.replace('PCR Adapter 384 Well and ', '')
                loc = 'PCR96WellAdapter_CoverSite_1'     # yes, it's correct
                pos = adapter_cnt['384 well']
            elif RackType.startswith('Alpaqua Magnum 96 Well and '):
                adapter_cnt['96 well magnet'] += 1
                # adding adapter to labware table
                df_labware.append({'labware_name' : 'Magnet for {}'.format(RackLabel),
                                   'labware_type' : 'Alpaqua Magnum FLX 96 well',
                                   'target_location' : loc,
                                   'target_position' : pos,
                                   'target_location_prompt' : loc,
                                   'target_position_prompt' : pos})
                # editing info for plate on adapter
                RackType = RackType.replace('Alpaqua Magnum 96 Well and ', '')
                loc = '96MicroplateSkirted_CoverSite_6'
                pos = adapter_cnt['96 well magnet'] 
            # if labware type is trough, ordering: 6-10,1-5
            try:
                lw_cat = v['category']
            except KeyError:
                raise KeyError('Cannot find "category" for labware: {}'.format(RackLabel))
            if lw_cat == 'trough':
                if pos <= 5:
                    pos += 5
                elif pos >= 6:
                    pos -= 5
                pos_prompt = pos            
            # creating table entry            
            df_labware.append({'labware_name' : RackLabel,
                               'labware_type' : RackType,
                               'target_location' : loc,
                               'target_position' : pos,
                               'target_location_prompt' : loc_prompt,
                               'target_position_prompt' : pos_prompt})  # txt = what is prompted for user
            
        # covert to dataframe & return
        df_tips = pd.DataFrame.from_dict(df_tips)
        df_labware = pd.DataFrame.from_dict(df_labware)
        ## ordering labware dataframe
        df_labware.sort_values(by=['target_location_prompt', 'target_position_prompt' ,'labware_type'],
                               inplace=True, ascending=[True, True, False])
        # return
        return pd.concat([df_tips, df_labware])
    
    def _next(self, labware, loc_tracker, keep_empty=True):
        """Getting nest target position for target location
        """
        # what the possible 
        psbl_targets = labware['target_location']
        # which of the possible target positions to use?
        target = None
        for x in psbl_targets:
            try:                
                target = self.target_position[x]
                target = [x, target]
                break
            except KeyError:
                continue            
        if target is None:
            return None
        # target variables
        try:
            position_count = target[1]['position_count']
        except KeyError:
            position_count = 1
        try:
            keep_empty_list = target[1]['keep_empty']
        except KeyError:
            keep_empty_list = None
        # setting target location
        while 1:
            # adding to position counter
            try:
                loc_tracker[target[0]] += 1
            except KeyError:
                loc_tracker[target[0]] = 1

            # checking if position is available
            if loc_tracker[target[0]] >= position_count:  
                msg = 'Not enough positions for target: "{}"'
                raise ValueError(msg.format(target[0]))
            elif (keep_empty == True and
                  loc_tracker[target[0]] in keep_empty_list):
                continue
            else:
                return target[0], loc_tracker[target[0]]                           
        return None

    def _add_labware(self, cmd, gwl):
        """Adding labware (no tip boxes) to self
        """
        if isinstance(cmd, Fluent.Reagent_distribution):
            assert cmd.SrcRackLabel is not None
            assert cmd.DestRackLabel is not None
            assert cmd.SrcRackType is not None
            assert cmd.DestRackType is not None
            # labware info
            self.labware[cmd.SrcRackLabel] = gwl.db.get_labware(cmd.SrcRackType)
            self.labware[cmd.DestRackLabel] = gwl.db.get_labware(cmd.DestRackType)
            # labware order in the gwl object
            try:
                _ = self.labware_order[cmd.SrcRackLabel]
            except KeyError:
                self.labware_order[cmd.SrcRackLabel] = len(self.labware_order.keys())
            try:
                _ = self.labware_order[cmd.DestRackLabel]
            except KeyError:
                self.labware_order[cmd.DestRackLabel] = len(self.labware_order.keys())            
        else:
            try:
                RackLabel = cmd.RackLabel
            except AttributeError:
                return None
            try:
                RackType = cmd.RackType
            except AttributeError:
                return None
            # labware info
            self.labware[RackLabel] = gwl.db.get_labware(RackType)
            # labware order in the gwl
            try:
                _ = self.labware_order[RackLabel]
            except KeyError:
                self.labware_order[RackLabel] = len(self.labware_order.keys()) 
    
    def _add_tip_boxes(self, gwl):
        """Adding tip boxes to self
        """
        # getting tip boxes from count
        for TipType,count in self.tip_count.items():
            tip_box = gwl.db.get_tip_box(TipType)
            d = gwl.db.get_labware(tip_box)
            try:
                wells = d['wells']
            except KeyError:
                msg = '"wells" key not found for labware: "{}"'
                raise KeyError(msg.format(tip_box))
            # number of tip boxes for the well
            n_boxes = int(round(count / wells + 0.5,0))
            for i in range(n_boxes):
                tip_box_label = '{0}[{1:0>3}]'.format(tip_box, i + 1)
                self.tip_boxes[tip_box_label] = [i, gwl.db.get_labware(tip_box)]
                        
    def _count_tips(self, commands):
        """Counting all tip usage in gwl commands and adding to self.
        Tips can be re-used if no waste command provided, so counting number 
        of Asp-Waste found.
        All Asp commands lacking a TipType will be skipped
        """
        TipType = None
        for cmd in commands:
            if isinstance(cmd, Fluent.Waste):
                # adding tip to count
                try:
                    self.tip_count[TipType] += 1
                except KeyError:
                    self.tip_count[TipType] = 1
            if isinstance(cmd, Fluent.Aspirate):
                # getting tip type for aspirate
                try:
                    TipType = cmd.TipType
                except AttributeError:
                    TipType = None
            if isinstance(cmd, Fluent.Reagent_distribution):
                # all tips used
                assert cmd.TipType is not None
                try: 
                    self.tip_count[cmd.TipType] += 8    # TODO: more precise
                except KeyError:
                    self.tip_count[cmd.TipType] = 8    # TODO: more precise
                
                    
class worktable_tracker():
    """Keeping track of available target positions for each target location

    NOTE: depreciated
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

    NOTE: depreciated
    """
    def __init__(self, tip_types={}):
        self.tip_types = self.tip_type_setter(tip_types)
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
        # adding tip to tip count; must have tip types
        if add_tip == True:
            tip_size = None
            for k,v in self.tip_types.items():
                try:
                    if x.Volume < k * 1.05 and v is not None:
                        tip_size = k
                except KeyError:
                    pass
            if tip_size is not None:
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

    def tip_for_volume(self, volume):
        """Return available tip type for volume (ul).
        Using smallest available tip type for volume.
        """
        for k in sorted(self.tip_types.keys()):
            if k > volume * 1.05 and self.tip_types[k] is not None:
                return self.tip_types[k]
        return ''            
        
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
        # return
        return df


    def tip_type_setter(self, value):
        for k,v in value.items():
            if v.lower() == 'none':
                value[k] = None
            else:
                try:
                    TIP_TYPE_DB[v]
                except KeyError:
                    msg = 'ERROR: tip type "{}" not recognized'
                    print(msg.format(v), file=sys.stderr)
                    sys.exit(1)
        self.tip_type = value
                        

# main
if __name__ == '__main__':
    pass
