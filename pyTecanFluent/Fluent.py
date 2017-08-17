# -*- coding: utf-8 -*-

# import
from __future__ import print_function
import os
import sys
import numpy as np


def xstr(x):
    if x is None:
        return ''
    else:
        return x

def _psbl_liq_cls():
    """Returns a set of possible liquid classes available for Fluent
    """
    x = ('Water Free Multi', 'Water Free Single',
         'MasterMix Free Multi', 'MasterMix Free Single', 
         'Ethanol Free Multi', 'Ethanol Free Single', 
         'DMSO Free Multi', 'DMSO Free Single', 
         'Serum Free Multi', 'Serum Free Single',
         'Water Contact Wet Multi', 'Water Contact Wet Single',
         'Water Mix')
    return x

class asp_disp():
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
    ForceRack
    MinDetected
    """

    def __init__(self):
        self._ID = ''
        # aspirate parameters
        self.RackLabel = None
        self.RackID = None
        self.RackType = None
        self._Position = 1
        self.TubeID = None
        self.Volume = None
        self._LiquidClass = 'Water Free Single'
        self.TipType = None
        self.TipMask = None
        self.key_order = ['_ID',
                          'RackLabel', 'RackID', 'RackType',
                          'Position', 'TubeID', 'Volume',
                          'LiquidClass', 'TipType', 'TipMask']
        self.psbl_liq_cls = _psbl_liq_cls()


    def cmd(self):
        # list of values in correct order
        vals = [getattr(self, x) for x in self.key_order]
        # None to blank string
        vals = [xstr(x) for x in vals]
        # convert all to strings
        vals = [str(x) for x in vals]
        # return
        return ';'.join(vals)

    def liquid_classes(self):
        x = '\n,'.join(list(self.psbl_liq_cls))
        print(x)        

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
        if value not in self.psbl_liq_cls:
            msg = 'Liquid class "{}" not allowed'
            raise TypeError(msg.format(value))
        self._LiquidClass = value



class aspirate(asp_disp):
    def __init__(self):
        asp_disp.__init__(self)
        self._ID = 'A'


class dispense(asp_disp):
    def __init__(self):
        asp_disp.__init__(self)
        self._ID = 'D'


class multi_disp():
    """Commands for aliquoting mastermix
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
        self._LiquidClass = 'Water Free Multi'
        self.NoOfMultiDisp = 2
        self.psbl_liq_cls = _psbl_liq_cls()

    def xstr(self, x):
        if x is None:
            return ''
        else:
            return x

    def cmd(self):
        # volume as interable
        if hasattr(self.Volume, '__iter__'):
            self.Volumes = self.Volume
        else:
            self.Volumes = [self.Volume] * len(self.DestPositions)
        # each multi-disp
        steps = []
        for i in range(0, self.SampleCount, self.NoOfMultiDisp):
            # number of dispenses
            if self.SampleCount - i < self.NoOfMultiDisp:
                n_disp = self.SampleCount - i
            else:
                n_disp = self.NoOfMultiDisp
            # single-asp
            asp = aspirate()
            asp.RackLabel = self.SrcRackLabel
            asp.Position = self.SrcPosition
            asp.Volume = sum(self.Volumes[i:(i+n_disp)])   
            asp.Volume = round(asp.Volume, 2)
            asp.LiquidClass = self.LiquidClass
            steps.append(asp.cmd())
            # multi-disp
            for ii in range(n_disp):
                disp = dispense()
                disp.RackLabel = self.DestRackLabel[i+ii]
                disp.Position = self.DestPositions[i+ii]
                disp.Volume = round(self.Volumes[i+ii], 2)
                disp.LiquidClass = self.LiquidClass
                steps.append(disp.cmd())
            steps.append('W;')
                
        return '\n'.join(steps)

    @property
    def LiquidClass(self):
        return self._LiquidClass

    @LiquidClass.setter
    def LiquidClass(self, value):
        if value not in self.psbl_liq_cls:
            msg = 'Liquid class "{}" not allowed'
            raise TypeError(msg.format(value))
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
        assert min(values) > 0, 'Min position is 1'
        self._DestPositions = values 

    @property
    def SampleCount(self):
        return len(self._DestPositions)



class reagent_distribution():
    """Commands for aliquoting mastermix
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

    def xstr(self, x):
        if x is None:
            return ''
        else:
            return x

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
