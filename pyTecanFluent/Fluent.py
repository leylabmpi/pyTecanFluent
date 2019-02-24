from __future__ import print_function

# import
import os
import sys
import re
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
    """Database of FluentControl labware, tip types, liquid classes, etc.
    Database files are stored in JSON format
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

    def RackTypes(self):
        return list(self.labware.keys())
            
    def get_labware(self, value):
        try:
            d = self.labware[value]
            d['RackType'] = value
            return d
        except KeyError:
            msg = 'Labware not in database: "{}"'
            raise KeyError(msg.format(value))

    def get_labware_wells(self, value):
        try:
            return self.labware[value]['wells']
        except KeyError:
            msg = 'Cannot find number of wells for labware: "{}"'
            raise KeyError(msg.format(value))                

    def get_labware_max_volume(self, value):
        try:
            return self.labware[value]['max_volume']
        except KeyError:
            msg = 'Cannot find max_volume for labware: "{}"'
            raise KeyError(msg.format(value))        

    def get_labware_allowed_tips(self, value):
        try:
            return self.labware[value]['allowed_tips']
        except KeyError:
            return None
        
    def get_tip_type(self, value):
        if value is None:
            return None
        try:
            return self.tip_type[value]
        except KeyError:
            msg = 'Tip type not in databse: "{}"'
            raise KeyError(msg.format(value))            

    def get_tip_volume(self, value):
        """Max volume of liquid that could fit in the tip if the liquid
        class doesn't require air volume"""
        try:
            return self.tip_type[value]['volume']
        except KeyError:
            msg = 'Cannot find volume for tip type: "{}"'
            raise KeyError(msg.format(value))            

    def get_tip_DTH_volume(self, value):
        """Max volume of liquid used for dynamic tip handling
        """
        try:
            return self.tip_type[value]['DTH']
        except KeyError:
            msg = 'Cannot find DTH (dynamic tip handling) volume for tip type: "{}"'
            raise KeyError(msg.format(value))           
        
    def get_tip_box(self, value):
        try:
            return self.tip_type[value]['tip_box']
        except KeyError:
            msg = 'Cannot find tip-box for tip type: "{}"'
            raise KeyError(msg.format(value))            

    def get_liquid_class(self, value):
        try:
            return self.liquid_class[value]
        except KeyError:
            msg = 'Liquid class not in database: "{}"'
            raise KeyError(msg.format(value))            
        
        
class gwl(object):
    """Class for storing gwl commands
    """
    def __init__(self, TipTypes=None):
        self.db = db()
        self.TipTypes = TipTypes
        self.last_asp = None
        self.commands = []

    def add(self, obj, default_liq_cls='Water Free Single'):
        """Adding gwl commands ('obj') to list of commands.
        A TipType will be added, but this is just used for counting tips later on
        """
        # assertions
        ## check values for asp/disp commands
        if isinstance(obj, Aspirate) or isinstance(obj, Dispense):
            assert obj.RackType is not None
            # check that RackType is in database
            labware = self.db.get_labware(obj.RackType)
            # check that tubes have target_position of 1 (only 1 position per tube)
            if 'eppendorf' in set(labware['target_location']) and int(obj.Position) != 1:
                obj.Position = 1
                
            # check that liquid class is in database (or use default)
            if self.LiquidClass_exists(obj.LiquidClass) is False:
                obj.LiquidClass = default_liq_cls
                
        # adding tip type for command based on volume; used for counting
        if isinstance(obj, Aspirate):
            obj.TipType = self.set_TipType(obj.Volume, obj.RackType)
            self.last_asp = obj
        elif isinstance(obj, Reagent_distribution):
            obj.TipType = self.set_TipType(obj.volume_per_aspirate())
        # dispense must have same liquid class as previous aspirate
        if isinstance(obj, Dispense):
            assert self.last_asp is not None
            obj.LiquidClass = self.last_asp.LiquidClass
                        
        # appending to list of commands
        self.commands.append(obj)

    def write(self, file_obj):
        """Writing out gwl file.
        Commands written in the order of addition
        """
        try:
            outF = open(file_obj, 'w')
        except IOError:
            outF = file_obj
            
        for x in self.commands:
            outF.write(x.cmd() + '\n')
        outF.close()
        
    def set_TipType(self, volume, racktype=None):
        """Setting which tip will be used.
        Tip selection based on dynamic tip handling (DTH) volume specified in the database
        """
        volume = float(volume)

        # check if racktype has min-volume tip
        ## if yes, set volume to that, which sets tip type
        if racktype is not None:
            tip_sizes = self.db.get_labware_allowed_tips(racktype)
            if tip_sizes is not None and min(tip_sizes) > volume:
                volume = min(tip_sizes) * 0.75    # WARNING: 0.75 is a hack!

        # setting tip type
        assert self.TipTypes is not None
        func = lambda x: (x[1]['DTH'],x[0])
        for k,v in sorted(self.TipTypes.items(), key=func):
            if v['DTH'] - volume > 0:
                return k            
        msg = 'No TipType DTH value greater than {}'
        raise ValueError(msg.format(volume))
            
    def TipType_exists(self, volume, warn=False):
        """Does TipType exist?
        """
        try:
            self.TipTypes[volume]
            return True
        except KeyError:
            if warn is True:
                msg = 'No TipType exists for volume: "{}"'
                print(msg.format(volume), file=sys.stderr)
            return False

    def get_DTH_volumes(self):
        """Getting the dynamic tip handling max volume for each TipType in gwl object
        """
        d = {}
        for TipType in self.TipTypes:
            d[TipType] = self.db.get_tip_DTH_volume(TipType)
        return d
        
    def RackType_count(self, RackType):
        """Counting labware with same RackType (different RackLabel, same RackType)
        """
        cnt = {}
        for x in self.commands:
            # getting possible target locations for RackType
            try:
                cmd_RackLabel = x.RackLabel
                cmd_RackType = x.RackType
            except AttributeError:
                continue
            if RackType == cmd_RackType:
                cnt[cmd_RackLabel] = 1
        return len(cnt.keys())
                
    def LiquidClass_exists(self, liquid_class, warn=False):
        """Check that liquid class exists in the database.
        """
        try:
            self.db.liquid_class[liquid_class]
            return True
        except KeyError:
            if warn is True:
                msg = 'Liquid class does not exist: "{}"'
                print(msg.format(liquid_class), file=sys.stderr)
            return False

    def list_commands(self):
        """List the commands added to the instance
        """
        return(self.commands)
        
    # get/set the available tip types (& their attributes) 
    @property
    def TipTypes(self):
        return self._TipTypes
    @TipTypes.setter
    def TipTypes(self, types):
        if types is None:
            self._TipTypes = None
        else:
            types = [x for x in types if x is not None]  
            self._TipTypes = {x:self.db.get_tip_type(x) for x in types}
                
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
        self.db = db()
        # aspirate parameters
        self.RackLabel = RackLabel
        self.RackID = RackID
        self.RackType = RackType
        self.Position = Position
        self.TubeID = TubeID
        self.Volume = Volume
        self.LiquidClass = LiquidClass
        self.TipType = TipType        # doesn't actually work!
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
        self.db.get_liquid_class(value)
        self._LiquidClass = value
        
class Aspirate(asp_disp):
    """gwl aspirate command: "A;"
    """
    def __init__(self):
        asp_disp.__init__(self)
        self._ID = 'A'

class Dispense(asp_disp):
    """gwl dispense command: "D;"
    """
    def __init__(self):
        asp_disp.__init__(self)
        self._ID = 'D'

class Comment():
    """gwl comment command: "C;"
    """
    def __init__(self, comment=''):
        self.comment = comment

    def cmd(self):
        return 'C;' + self.comment.lstrip('C;')

class Waste():
    """gwl waste command: "W;"
    Used for ejecting tip. 
    """
    def __init__(self):
        pass

    def cmd(self):
        return 'W;'

class Flush():
    """gwl flush command: "F;"
    This is useful for re-using tips.
    The flush command will flush the extra
    volume back into the source labware.
    Note: the source labware must have the variable: "IsFCAWaste = True"
    Using the following worklist command structure : "A; D; F; ... W; B;"
    """
    def __init__(self):
        pass

    def cmd(self):
        return 'F;'
    
class Break():
    """gwl break command: "B;"
    Useful for forcing tip ejects and preventing errors 
    transitioning between pipetting tasks (eg., mastermix & water)
    """
    def __init__(self):
        pass

    def cmd(self):
        return 'B;'
    
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
        self.LiquidClass = 'Water Free Multi'
        self.NoOfMultiDisp = 2

    def add(self, gwl, disp_frac):
        # volume as iterable
        if hasattr(self.Volume, '__iter__'):
            self.Volumes = self.Volume
        else:
            self.Volumes = [self.Volume] * len(self.DestPositions)
                    
        # each multi-disp
        sample_cnt = 0
        while 1:
            # single-asp
            asp = Aspirate()
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
                disp = Dispense()
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
            asp.Volume = round(sum([x.Volume for x in dispenses]) * (1-disp_frac+1), 2)

            # appending to gwl-obj
            gwl.add(asp)
            for x in dispenses:
                gwl.add(x)
            gwl.add(Waste())                       

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

    
class Reagent_distribution(object):
    """
    NOTE: Using this command requires a labware designated for FCA waste
    (eg., a trough or tube)
    The labware must have "IsFcaLiquidWaste" == True in the "Custom Attributes"

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
    Direction = Which way to pipette? Default:0 (column-wise)
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
        self.DestPosStart = 1
        self.DestPosEnd = 1
        # other
        self.Volume = 1
        self.LiquidClass = 'Water Free Multi'
        self.NoOfDiTiReuses = 1
        self.NoOfMultiDisp = 5
        self.Direction = 0
        self.ExcludedDestWell = None
        self.TipType = None      # just used for counting tips
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

    def volume_per_aspirate(self):
        """Get the volume per aspiration for the multi-dispense
        """
        return self.Volume * self.NoOfMultiDisp

# main
if __name__ == '__main__':
    pass
