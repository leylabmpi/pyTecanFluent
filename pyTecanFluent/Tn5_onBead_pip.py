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
def batch(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]

def pip_mastermix(df_map, gwl,  mm_labware_type='25ml_1 waste',
                  mm_volume=13.1, liq_cls='MasterMix Free Single',
                  sup_volume=100, sup_rm_liq_cls='Water Free Single'):
    """
    Writing worklist commands for aliquoting mastermix.
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
#    x = cycle(range(8))
#    df['CHANNEL_ORDER'] = [next(x) for y in range(df.shape[0])]
#    x = cycle(range(n_tip_reuse))
#    df['TIP_BATCH'] = Utils.tip_batch(df['CHANNEL_ORDER'], n_tip_reuse)
#    df.sort_values(by=['TIP_BATCH',
#                       'CHANNEL_ORDER',
#                       'TECAN_dest_target_position'], inplace=True)
#    df.reset_index(inplace=True)
#    print(df); sys.exit()

    # determining sup. asp-disp cycles
    n_cycles = 1 #int(round(sup_volume / 192 + 0.49999,0))
    sup_volume = round(sup_volume / n_cycles,1)

    # for each Sample-PCR, write out asp/dispense commands
    for x in batch(range(df.shape[0]), 8):
        for i in range(n_cycles):
            for ii in x:
                # removing supernatant
                ## aspiration
                asp = Fluent.Aspirate()
                asp.RackLabel = df.loc[ii,'TECAN_dest_labware_name']
                asp.RackType = df.loc[ii,'TECAN_dest_labware_type']
                asp.Position = df.loc[ii,'TECAN_dest_target_position']
                asp.Volume = sup_volume
                asp.LiquidClass = sup_rm_liq_cls
                gwl.add(asp)
                ## dispensing
                disp = Fluent.Dispense()
                disp.RackLabel = '100ml_waste[001]'
                disp.RackType = '100ml_1'
                disp.Position = 1
                disp.Volume = sup_volume
                disp.LiquidClass = sup_rm_liq_cls
                gwl.add(disp)
                if i + 1 == n_cycles or len([y for y in x]) < 8:
                    gwl.add(Fluent.Waste())
                else:
                    gwl.add(Fluent.Flush())
        gwl.add(Fluent.Break())
        for i in x:    
            # aliquoting mastermix
            ## aspiration
            asp = Fluent.Aspirate()
            asp.RackLabel = 'Mastermix[{0:0>3}]'.format(1)
            asp.RackType = mm_labware_type
            asp.Position = 1
            asp.Volume = mm_volume
            asp.LiquidClass = liq_cls
            gwl.add(asp)
            ## dispensing
            disp = Fluent.Dispense()
            disp.RackLabel = df.loc[i,'TECAN_dest_labware_name']
            disp.RackType = df.loc[i,'TECAN_dest_labware_type']
            disp.Position = df.loc[i,'TECAN_dest_target_position']
            disp.Volume = mm_volume
            disp.LiquidClass = liq_cls
            gwl.add(disp)
            ## tip waste
            gwl.add(Fluent.Waste())
    gwl.add(Fluent.Break())
                
def pip_primer(i, gwl, df_map, primer_labware_name, 
               primer_labware_type, primer_target_position,
               primer_volume, liq_cls):
    """
    Pipetting a single primer
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
    """
    Commands for aliquoting primers
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
