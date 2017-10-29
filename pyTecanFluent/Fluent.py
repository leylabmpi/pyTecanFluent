from __future__ import print_function

# import
import os
import sys
import json
import collections
import pkg_resources
import numpy as np
import pandas as pd


#-- notes on gwl file format --#
# A;RackLabel;RackID;RackType;Position;TubeID;Volume;LiquidClass;Tip Type;TipMask;ForcedRackType
## 11 fields; 10 semicolons
# D;RackLabel;RackID;RackType;Position;TubeID;Volume;LiquidClass;Tip Type;TipMask;ForcedRackType
## 11 fields; 10 semicolons

def xstr(x):
    if x is None:
        return ''
    else:
        return x


class db(object):
    """
    Database of FluentControl labware, tip types, liquid classes, etc
    """
    def __init__(self):
        d = os.path.split(__file__)[0]
        self.database_dir = os.path.join(d, 'database')
        # labware
        f = os.path.join(self.database_dir, 'labware.json')
        with open(f) as inF:
            self.labware = json.load(inF)
        # tip type
        f = os.path.join(self.database_dir, 'tip_type.json')
        with open(f) as inF:
            self.tip_type = json.load(inF)
        # liquid class
        f = os.path.join(self.database_dir, 'liquid_class.json')
        with open(f) as inF:
            self.liquid_class = json.load(inF)

    def labware_keys(self):
        return self.labware.keys()
    
    def get_labware(self, value):
        try:
            return self.labware[value]
        except KeyError:
            msg = 'Cannot find labware "{}"'
            raise KeyError(msg.format(value))

    def get_labware_max_volume(self, value):
        try:
            return self.labware[value]['max_volume']
        except KeyError:
            msg = 'Cannot find max_volume for labware "{}"'
            raise KeyError(msg.format(value))        
        
    def get_tip_type(self, value):
        try:
            return self.tip_type[value]
        except KeyError:
            msg = 'Cannot find tip type "{}"'
            raise KeyError(msg.format(value))            

    def get_tip_volume(self, value):
        try:
            return self.tip_type[value]['volume']
        except KeyError:
            msg = 'Cannot find volume for tip type "{}"'
            raise KeyError(msg.format(value))            

    def get_tip_box(self, value):
        try:
            return self.tip_type[value]['tip_box']
        except KeyError:
            msg = 'Cannot find tip-box for tip type "{}"'
            raise KeyError(msg.format(value))            

    def get_liquid_class(self, value):
        try:
            return self.liquid_class[value]
        except KeyError:
            msg = 'Cannot find liquid class "{}"'
            raise KeyError(msg.format(value))            
        

class labware(object):
    """Class for summarizing labware in a gwl object
    """
    def __init__(self):
        self.tip_count = {}
        self.tip_boxes = {}
        self.labware = collections.OrderedDict()
        # target position
        d = os.path.join(os.path.split(__file__)[0],  'database')
        f = os.path.join(d, 'target_position.json')
        with open(f) as inF:
            self.target_position = json.load(inF)
                
    def add_gwl(self, gwl):
        """Adding labware from gwl object to labware object.
        Note: this can be used to sum up labware from multiple gwl objects.
        """
        for cmd in gwl.commands:
            # counting tips
            self._count_tips(cmd)
            # labware IDs
            self._add_labware(cmd, gwl)
        # summing up tip boxes        
        self._add_tip_boxes(gwl)

    def table(self):
        """Creating pandas dataframe of labware
        columns: labware_name, labware_type,target_location,target_position
        """
        #print(self.tip_boxes)
        #print(self.labware)
        #print(self.target_position)

        # init target position counters
        loc_tracker = {}
        
        # init dataframe
        cols = ['labware_name', 'labware_type',
                'target_location', 'target_position']
        df = pd.DataFrame(columns=cols)
        
        # adding tip boxes
        #for k,v in sorted(self.tip_boxes.iteritems(),
        #                  key=lambda (k,v): (v,k),
        #                  reverse=True):
        for k,v in sorted(self.tip_boxes.items()):
            loc = self._next_location(v, loc_tracker)
            if loc is None:
                msg = 'No possible target location for labware: "{}"'
                raise ValueError(msg.format(k))
            print(loc)

    def _next_location(self, labware, loc_tracker):
        # what the possible 
        psbl_targets = labware['target_location']
        # which of the possible target positions to use?
        target = None
        for x in psbl_targets:
            try:                
                target = self.target_position[x]
                target = [x, target]
            except KeyError:
                continue
        if target is None:
            return None
        # setting target location
        ## TODO: use while loop to account for boards and running out of positions
        try:
            loc_tracker[target[0]] += 1
        except KeyError:
            loc_tracker[target[0]] = 1
        
        return loc_tracker[target[0]]
                    
        
    def _add_tip_boxes(self, gwl):
        """Adding tip boxes to labware
        """
        # getting tip boxes from count
        for k in self.tip_count.keys():
            tip_box = gwl._db.get_tip_box(k)
            self.tip_boxes[tip_box] = gwl._db.get_labware(tip_box)

    def _add_labware(self, cmd, gwl):
        try:
            RackLabel = cmd.RackLabel
        except AttributeError:
            return None
        try:
            RackType = cmd.RackType
        except AttributeError:
            return None
        self.labware[RackLabel] = gwl._db.get_labware(RackType)
            
    def _count_tips(self, cmd):
        """Counting all tips in gwl commands
        """
        # Just count aspirations
        if not isinstance(cmd, aspirate):
            return None
        # Tip type count
        try:
            TipType = cmd.TipType
        except AttributeError:
            return None
        try:
            self.tip_count[TipType] += 1
        except KeyError:
            self.tip_count[TipType] = 1

            
class gwl(object):
    """Class for storing gwl commands
    """
    def __init__(self, TipTypes=None):
        self._db = db()
        self.TipTypes = TipTypes
        self._last_asp_TipType = None
        self.commands = []

    def add(self, obj, force_tip=True):
        """Adding gwl commands ('obj') to list of commands
        """
        # assertions
        if isinstance(obj, aspirate) or isinstance(obj, dispense):
            assert obj.RackType is not None
        # forcing usage of particular tip type
        if force_tip is True:
            if isinstance(obj, dispense):
                assert self._last_asp_TipType is not None
                obj.TipType = self._last_asp_TipType
            elif isinstance(obj, aspirate):
                obj.TipType = self.set_TipType(obj.Volume)
                self._last_asp_TipType = obj.TipType
        # appending to list of commands
        self.commands.append(obj)

    def write(self, file_obj):
        """Writing out gwl fileC
        Commands written in the order of addition
        """
        try:
            outF = open(file_obj, 'w')
        except IOError:
            outF = file_obj
        
        for x in self.commands:
            outF.write(x.cmd() + '\n')
        outF.close()
                    
    def set_TipType(self, volume):
        """Setting which tip is used for the command
        """
        if self.TipTypes is None:
            return None
        for k in sorted(self.TipTypes.keys()):
            if k > volume * 1.05 and self.TipTypes[k] is not None:
                return self.TipTypes[k]
        return None
    
    @property
    def TipTypes(self):
        return self._TipTypes
    @TipTypes.setter
    def TipTypes(self, value):
        for k,v in value.items():
            # key should be integer or float
            try:
                k = float(k)
            except ValueError:
                msg = 'TypeType key "{}" cannot be converted to float'
                raise ValueError(msg.format(k))
            # values should be valid tip type            
            self._db.get_tip_type(v)                
        self._TipTypes = value

                
class asp_disp(object):
    """Commands for aliquoting mastermix
    *Parameters*
    RackLabel
    RackID
    RackType
    Position
    TubeID
    Volume
    LiquidClass 
    TipType
    TipMask
    ForceRackType
    """
    def __init__(self, RackLabel=None, RackID=None, RackType=None,
                 Position=1, TubeID=None, Volume=None,
                 LiquidClass = 'Water Free Single', TipType=None,
                 TipMask=None, ForceRack=None):
        self._ID = ''
        self._db = db()
        # aspirate parameters
        self.RackLabel = RackLabel
        self.RackID = RackID
        self.RackType = RackType
        self.Position = Position
        self.TubeID = TubeID
        self.Volume = Volume
        self.LiquidClass = LiquidClass
        self.TipType = TipType
        self.TipMask = TipMask
        self.ForceRackType = ForceRack
        self.field_order = ['_ID',
                            'RackLabel', 'RackID', 'RackType',
                            'Position', 'TubeID', 'Volume',
                            'LiquidClass', 'TipType', 'TipMask',
                            'ForceRackType']
        
    def cmd(self):
        # assertions
        assert self.RackLabel is not None, 'RackLabel cannot be None'
        assert self.RackType is not None, 'RackType cannot be None'
        # list of values in correct order
        vals = [getattr(self, x) for x in self.field_order]
        # None to blank string
        vals = [xstr(x) for x in vals]
        # convert all to strings
        vals = [str(x) for x in vals]
        # return
        return ';'.join(vals)

    @property
    def Position(self):
        return self._Position

    @Position.setter
    def Position(self, value):
        self._Position = int(value)

    @property
    def LiquidClass(self):
        return self._LiquidClass

    @LiquidClass.setter
    def LiquidClass(self, value):
        self._db.get_liquid_class(value)
        self._LiquidClass = value
        
class aspirate(asp_disp):
    """gwl aspirate command: "A;"
    """
    def __init__(self):
        asp_disp.__init__(self)
        self._ID = 'A'

class dispense(asp_disp):
    """gwl dispense command: "D;"
    """
    def __init__(self):
        asp_disp.__init__(self)
        self._ID = 'D'

class comment():
    """gwl comment command: "C;"
    """
    def __init__(self, comment=''):
        self.comment = comment

    def cmd(self):
        return 'C;' + self.comment.lstrip('C;')

class waste():
    """gwl waste command: "W;"
    """
    def __init__(self):
        pass

    def cmd(self):
        return 'W;'
    
class multi_disp(object):
    """Commands for aliquoting reagent to multiple labware positions
    *AspirateParameters*
    SrcRackLabel
    SrcRackID
    SrcRackType
    SrcPosition
    *DispenseParameters*
    DestRackLabel
    DestRackID
    DestRackType
    DestPosStart
    DestPosEnd
    *Samples*
    SampleCount
    *Other*
    Volume = How much volume per dispense?
    LiquidClass = Which liquid class to use? Default: 'Water Free Multi'
    NoOfMultiDisp = How many multi-dispenses?
    Labware_tracker = Labware object that tracks what labware is needed
    Returns
    * string of commands
    """

    def __init__(self):
        self._ID = 'R'
        # aspirate parameters
        self.SrcRackLabel = None
        self.SrcRackID = None
        self.SrcRackType = None
        self.SrcPosition = 1
        # dispense parameters
        self.DestRackLabel = []
        self.DestRackID = []
        self.DestRackType = []
        self._DestPositions = [1]
        # other
        self.Volume = 1.0
        self.TipType = None
        self._LiquidClass = 'Water Free Multi'
        self.NoOfMultiDisp = 2
        #self.psbl_liq_cls = _psbl_liq_cls()

    def add(self, gwl):
        # volume as iterable
        if hasattr(self.Volume, '__iter__'):
            self.Volumes = self.Volume
        else:
            self.Volumes = [self.Volume] * len(self.DestPositions)
                    
        # each multi-disp
        sample_cnt = 0
        while 1:
            # single-asp
            asp = aspirate()
            asp.RackLabel = self.SrcRackLabel
            asp.RackType = self.SrcRackType
            asp.Position = self.SrcPosition
            asp.LiquidClass = self.LiquidClass
            # determining total volume for this asp            
            dispenses_tmp = 0
            sample_cnt_tmp = sample_cnt
            total_asp_volume = 0
            while 1:
                sample_cnt_tmp += 1
                # Total samples reached
                if sample_cnt_tmp > len(self.DestPositions):
                    break
                # Number of multi-disp reached
                if dispenses_tmp >= self.NoOfMultiDisp:
                    sample_cnt_tmp -= 1 
                    break
                # Skipping 0-volumes
                if self.Volumes[sample_cnt_tmp-1] <= 0:
                    continue
                disp_volume = round(self.Volumes[sample_cnt_tmp-1], 2)
                if disp_volume > 0:
                    total_asp_volume += disp_volume
                    dispenses_tmp += 1
            # loading dispenses
            dispenses = []
            while 1:
                sample_cnt += 1
                # Total samples reached
                if sample_cnt > len(self.DestPositions):
                    break
                # Number of multi-disp reached
                if len(dispenses) >= self.NoOfMultiDisp:
                    sample_cnt -= 1 
                    break
                # Skipping 0-volumes
                if self.Volumes[sample_cnt-1] <= 0:
                    continue 
                disp = dispense()
                disp.RackLabel = self.DestRackLabel[sample_cnt-1]
                disp.RackType = self.DestRackType[sample_cnt-1]
                disp.Position = self.DestPositions[sample_cnt-1]
                disp.Volume = round(self.Volumes[sample_cnt-1], 2)
                disp.LiquidClass = self.LiquidClass
                if disp.Volume > 0:
                    dispenses.append(disp)
            # break if no more dispenses
            if len(dispenses) <= 0:
                break
            # adding asp-disp cycle
            asp.Volume = round(sum([x.Volume for x in dispenses]) * 1.05, 1)

            # appending to gwl-obj
            gwl.add(asp)
            for x in dispenses:
                gwl.add(x)
            gwl.add(waste())                       

    @property
    def LiquidClass(self):
        return self._LiquidClass

    @LiquidClass.setter
    def LiquidClass(self, value):
        self._db.get_liquid_class(value)
        self._LiquidClass = value

    @property
    def DestPositions(self):
        return self._DestPositions

    @DestPositions.setter
    def DestPositions(self, values):
        try:
            values = values.tolist()
        except AttributeError:
            pass
        values = [int(x) for x in values]
        assert min(values) > 0, 'Min position is <= 0'
        self._DestPositions = values 

    @property
    def SampleCount(self):
        return len(self._DestPositions)

    
class reagent_distribution(object):
    """
    # NOTE: depreciated 

    Commands for aliquoting mastermix
    *AspirateParameters*
    SrcRackLabel
    SrcRackID
    SrcRackType
    SrcPosStart
    SrcPosEnd
    *DispenseParameters*
    DestRackLabel
    DestRackID
    DestRackType
    DestPosStart
    DestPosEnd
    *Other*
    Volume = How much volume to asp/disp?
    LiquidClass = Which liquid class to use? Default: 'Water Free Multi'
    NoOfDiTiReuses = How many times to reuse tips?
    NoOfMultiDisp = How many multi-dispenses?
    Direction = Which way to pipette? Default:0
    ExcludedDestWell = Semi-colon separated string of locations to exclude

    *WashParameters*
    None?

    # Example: R;100ml_2;;Trough 100ml;1;1;96 Well Skirted PCR[003];;96 Well Skirted PCR;1;96;20;Water Free Multi;1;5;0
    """

    def __init__(self, ):
        self._ID = 'R'
        # aspirate parameters
        self.SrcRackLabel = None
        self.SrcRackID = None
        self.SrcRackType = None
        self.SrcPosStart = 1
        self.SrcPosEnd = 1
        # dispense parameters
        self.DestRackLabel = None
        self.DestRackID = None
        self.DestRackType = None
        self.DestPositions = 1
        self.DestPosEnd = 1
        # other
        self.Volume = 1
        self.LiquidClass = 'Water Free Multi'
        self.NoOfDiTiReuses = 1
        self.NoOfMultiDisp = 5
        self.Direction = 0
        self.ExcludedDestWell = None
        self.key_order = ['_ID',
                          'SrcRackLabel', 'SrcRackID', 'SrcRackType',
                          'SrcPosStart', 'SrcPosEnd',
                          'DestRackLabel', 'DestRackID', 'DestRackType',
                          'DestPosStart', 'DestPosEnd',
                          'Volume', 'LiquidClass', 'NoOfDiTiReuses',
                          'NoOfMultiDisp', 'Direction', 'ExcludedDestWell']

    def cmd(self):
        # list of values in correct order
        vals = [getattr(self, x) for x in self.key_order]
        # None to blank string
        vals = [xstr(x) for x in vals]
        # convert all to strings
        vals = [str(x) for x in vals]
        # return
        return ';'.join(vals)



# main
if __name__ == '__main__':
    pass
