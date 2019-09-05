from __future__ import print_function
# import
## batteries
from itertools import product,cycle
## 3rd party
import numpy as np
## package
from pyTecanFluent import Utils
from pyTecanFluent import Fluent
from pyTecanFluent import Labware


# functions
def pip_Tn5_mastermix(df_map, gwl, mm_volume,
                      src_labware_type='2ml Eppendorf waste',
                      liq_cls='Water Free Single', n_tip_reuse=1):
    """Commands for aliquoting Tn5 mastermix (Tn5 + buffer + water)
    """    
    if mm_volume <= 0:
        return None

    # reordering to account for n_tip_reuse
    df_map = Utils.reorder_plate_n_tip_reuse(df_map, n_tip_reuse)

    # for each Sample, write out asp/dispense commands
    gwl.add(Fluent.Comment('Tn5 mastermix (Tn5 + buffer + water)'))
    for i in range(df_map.shape[0]):        
        # aspiration
        asp = Fluent.Aspirate()
        asp.RackLabel = 'Tn5_mastermix'
        asp.RackType = src_labware_type
        asp.Position = 1
        asp.Volume = mm_volume
        asp.LiquidClass = liq_cls
        gwl.add(asp)

        # dispensing
        disp = Fluent.Dispense()
        disp.RackLabel = df_map.loc[i,'TECAN_dest_labware_name']
        disp.RackType = df_map.loc[i,'TECAN_dest_labware_type']
        disp.Position = df_map.loc[i,'TECAN_dest_target_position']
        disp.Volume = mm_volume
        disp.LiquidClass = liq_cls
        gwl.add(disp)

        # waste
        if (i + 1) % n_tip_reuse == 0 or i + 1 == df_map.shape[0]:
            gwl.add(Fluent.Waste())
        
    # adding break
    gwl.add(Fluent.Break())

def pip_samples(df_map,  gwl, DNA_volume=1, liq_cls='Water Free Single', n_tip_reuse=1):
    """Commands for aliquoting samples to each PCR rxn
    """
    gwl.add(Fluent.Comment('Samples'))
    # for each Sample-PCR, write out asp/dispense commands
    for i in range(df_map.shape[0]):
        
        # aspiration
        asp = Fluent.Aspirate()
        asp.RackLabel = df_map.loc[i,'TECAN_sample_labware_name']
        asp.RackType = df_map.loc[i,'TECAN_sample_labware_type']
        asp.Position = df_map.loc[i,'TECAN_sample_target_position']
        asp.Volume = DNA_volume
        asp.LiquidClass = liq_cls
        gwl.add(asp)

        # dispensing
        disp = Fluent.Dispense()
        disp.RackLabel = df_map.loc[i,'TECAN_dest_labware_name']
        disp.RackType = df_map.loc[i,'TECAN_dest_labware_type']
        disp.Position = df_map.loc[i,'TECAN_dest_target_position']
        disp.Volume = DNA_volume
        disp.LiquidClass = liq_cls
        gwl.add(disp)

        # waste
        gwl.add(Fluent.Waste())
        
    # adding break
    gwl.add(Fluent.Break())

def pip_mastermix(df_map, gwl,  mm_labware_type='25ml_1 waste',
                  mm_volume=13.1, n_tip_reuse=6,
                  liq_cls='MasterMix Free Single'):
    """Writing worklist commands for aliquoting mastermix.
    Using 1-asp-multi-disp with 200 ul tips.
    Method:
    * calc max multi-dispense for 200 ul tips & mm_volume
    * for 1:n_dispense
      * determine how many disp left for channel (every 8th)
      * if n_disp_left < n_disp: n_disp = n_disp_left
      * calc total volume: n_disp * mm_volume
    """
    if mm_volume <= 0:
        return None    
    gwl.add(Fluent.Comment('MasterMix'))

    # copying df
    df = df_map.copy()
    x = cycle(range(8))
    df['CHANNEL_ORDER'] = [next(x) for y in range(df.shape[0])]
    x = cycle(range(n_tip_reuse))
    df['TIP_BATCH'] = Utils.tip_batch(df['CHANNEL_ORDER'], n_tip_reuse)
    df.sort_values(by=['TIP_BATCH',
                       'CHANNEL_ORDER',
                       'TECAN_dest_target_position'], inplace=True)
    df.reset_index(inplace=True)
    
    # for each Sample-PCR, write out asp/dispense commands
    for i in range(df.shape[0]):
        # aspiration
        asp = Fluent.Aspirate()
        asp.RackLabel = 'Mastermix[{0:0>3}]'.format(1)
        asp.RackType = mm_labware_type
        asp.Position = 1
        asp.Volume = mm_volume
        asp.LiquidClass = liq_cls
        gwl.add(asp)

        # dispensing
        disp = Fluent.Dispense()
        disp.RackLabel = df.loc[i,'TECAN_dest_labware_name']
        disp.RackType = df.loc[i,'TECAN_dest_labware_type']
        disp.Position = df.loc[i,'TECAN_dest_target_position']
        disp.Volume = mm_volume
        disp.LiquidClass = liq_cls
        gwl.add(disp)

        # tip waste
        if (i + 1) % n_tip_reuse == 0 or i + 1 == df.shape[0]:
            gwl.add(Fluent.Waste())
            
    # adding break
    gwl.add(Fluent.Break())

def pip_primer(i, gwl, df_map, primer_labware_name, 
               primer_labware_type, primer_target_position,
               primer_volume, liq_cls):
    """Pipetting a single primer
    """
    # aspiration
    asp = Fluent.Aspirate()
    asp.RackLabel = df_map.loc[i, primer_labware_name]
    asp.RackType = df_map.loc[i, primer_labware_type]
    asp.Position = df_map.loc[i, primer_target_position]
    asp.Volume = primer_volume
    asp.LiquidClass = liq_cls
    gwl.add(asp)

    # dispensing
    disp = Fluent.Dispense()
    disp.RackLabel = df_map.loc[i,'TECAN_dest_labware_name']
    disp.RackType = df_map.loc[i,'TECAN_dest_labware_type']
    disp.Position = df_map.loc[i,'TECAN_dest_target_position']
    disp.Volume = primer_volume
    asp.LiquidClass = liq_cls
    gwl.add(disp)

    # waste
    gwl.add(Fluent.Waste())
    
def pip_primers(df_map, gwl, prm_volume=0, liq_cls='Water Free Single'):
    """Commands for aliquoting primers
    """
    gwl.add(Fluent.Comment('Primers'))    
    for i in range(df_map.shape[0]):
        if prm_volume <= 0:
            continue
        pip_primer(i, gwl, df_map,
                   'TECAN_primer_labware_name', 
                   'TECAN_primer_labware_type',
                   'TECAN_primer_target_position',
                   prm_volume, liq_cls)
        
    # adding break
    gwl.add(Fluent.Break())
