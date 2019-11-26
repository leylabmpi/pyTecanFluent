# -*- coding: utf-8 -*-
# import
## batteries
from __future__ import print_function
import os
import sys
import argparse
from itertools import product,cycle
import string
## 3rd party
import numpy as np
import pandas as pd
## package
from pyTecanFluent import Utils
from pyTecanFluent import Fluent
from pyTecanFluent import Labware


# functions
def get_desc():
    desc = 'Create robot commands for qPCR assay setup'
    return desc

def parse_args(test_args=None, subparsers=None):
    # desc
    desc = get_desc()
    epi = """DESCRIPTION:
    Create a worklist file for the TECAN Fluent robot for qPCR setup.
    The input is an exported plate layout from the BioRad qPCR software.
    The file format should be Excel or CSV.
    Just create a plate layout for your experimet, then export and add some needed info:
      The following columns should also be added to the table:
      * "Sample labware name" 
         * labware name containing the sample (any name that you want)
         * Exmaple: "source plate"
      * "Sample labware type" 
         * labware type (must EXACTLY match an existing labware type)
         * Example: "96 Eppendorf TwinTec PCR"
      * "Sample location"  
         * location of sample in the source plate. 
         * numeric; column-wise indexing
      * "Sample volume"  
         * numeric; sample volume in ul
      * "MM name"  
         * Name of master mix for that sample
         * This allows for multiple master mixes per assay 
      * "MM volume"  
         * Volume of master mix in PCR rxn (ul)
      * "Water volume"  
         * Volume of water in PCR rxn (ul)
    
    Notes:
    * Sample locations in plates numbered are column-wise (left-to-right)
    * The setup file (input table) MUST have a header (capitalization doesn't matter)
    * All volumes are in ul.
    """
    if subparsers:
        parser = subparsers.add_parser('qPCR', description=desc, epilog=epi,
                                       formatter_class=argparse.RawTextHelpFormatter)
    else:
        parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                         formatter_class=argparse.RawTextHelpFormatter)

    # args
    ## I/O
    groupIO = parser.add_argument_group('I/O')
    groupIO.add_argument('setup', metavar='SetupFile', type=str,
                         help='An Excel or CSV file with experimental setup')
    groupIO.add_argument('--prefix', type=str, default='TECAN_qPCR',
                         help='Output file name prefix (default: %(default)s)')
    groupIO.add_argument('--format', type=str, default=None,
                         choices=[None, 'excel', 'csv', 'tsv'],
                         help='File format (excel, csv, or tsv). If not provided, the format is determined from the file extension (default: %(default)s)') 

    ## Source labware 
    src = parser.add_argument_group('Source labware')
    src.add_argument('--mm-type', type=str, default='2ml Eppendorf waste',
                      help='Mastermix labware type (default: %(default)s)')
    src.add_argument('--water-type', type=str, default='100ml_1 waste',
                      help='Water labware type (default: %(default)s)')
    
    ## Destination labware
    dest = parser.add_argument_group('Destination labware')
    dest.add_argument('--dest', type=str, default='Destination plate',
                      help='Destination plate labware name (default: %(default)s)')
    dest.add_argument('--dest-type', type=str, default='384 Well Biorad PCR',
                      help='Destination plate labware type (default: %(default)s)')

    # Liquid classes
    liq = parser.add_argument_group('Liquid classes')
    liq.add_argument('--mm-liq', type=str, default='MasterMix Free Multi Wall Disp',
                      help='Mastermix liquid class (default: %(default)s)')
    liq.add_argument('--samp-liq', type=str, default='Water Free Single Wall Disp',
                      help='Sample liquid class (default: %(default)s)')
    liq.add_argument('--water-liq', type=str, default='Water Free Single Wall Disp',
                      help='Water liquid class (default: %(default)s)')
    liq.add_argument('--n-tip-reuse', type=int, default=4,
                     help='Number of tip reuses for applicable reagents (default: %(default)s)')

    # Parse & return
    if test_args:
        args = parser.parse_args(test_args)
        return args

    return parser


def main(args=None):
    # Input
    if args is None:
        args = parse_args()
    check_args(args)

    # Load input assay setup table
    df_setup = load_setup(args.setup, 
                          file_format=args.format)
    
    # gwl object init 
    TipTypes = ['FCA, 1000ul SBS', 'FCA, 200ul SBS',
                'FCA, 50ul SBS', 'FCA, 10ul SBS']    
    gwl = Fluent.gwl(TipTypes)
    
    # adding sample/reagent destinations to setup table
    n_wells = gwl.db.get_labware_wells(args.dest_type)
    add_dest(df_setup, args.dest, args.dest_type, n_wells)
    df_setup = check_rack_labels(df_setup)
    
    # Reordering dest for optimal pipetting
    if n_wells == 384:
        df_setup = Utils.reorder_384well(df_setup, gwl,
                                       labware_name_col='dest_labware_name',
                                       labware_type_col='dest_labware_type',
                                       position_col='dest_target_position')
    elif n_wells == 96:
        df_setup.sort(['dest_target_position'], inplace=True)
    else:
        msg = 'Labware type "{}" not recognized'
        raise ValueError(msg.format(args.desttype))
    
    # Adding commands to gwl object
    pip_mastermixes(df_setup, gwl=gwl, 
                    src_labware_type=args.mm_type,
                    liq_cls=args.mm_liq,
                    n_tip_reuse=args.n_tip_reuse)
    
    ## Samples
    pip_samples(df_setup, gwl=gwl,
                liq_cls=args.samp_liq)
    
    ## Water
    pip_water(df_setup, gwl=gwl,
              src_labware_type=args.water_type,
              liq_cls=args.water_liq)

    ## writing out worklist (gwl) file
    gwl_file = args.prefix + '.gwl'
    gwl.write(gwl_file)
    
    # making labware table
    lw = Labware.labware()
    lw.add_gwl(gwl)
    lw_df = lw.table()
    lw_file = args.prefix + '_labware.txt'
    lw_df.to_csv(lw_file, sep='\t', index=False)

    # Creating report file
    report_file = args.prefix + '_report.txt'
    with open(report_file, 'w') as repFH:        
        MM_names = np.unique(df_setup['mm name'])
        for i,MM_name in enumerate(MM_names):
            df = df_setup.loc[df_setup['mm name'] == MM_name]        
            df.reset_index(inplace=True)
            write_report(df, MM_name=MM_name, outFH=repFH)
    
    # status on files written
    Utils.file_written(gwl_file)
    Utils.file_written(lw_file)
    Utils.file_written(report_file)
    
def check_args(args):
    """Checking user input
    """
    # special characters for namings
    args.dest = Utils.rm_special_chars(args.dest)
    
def load_setup(input_file, file_format=None, header=0):
    """Loading setup file (Excel, csv, or tab-delim)
    """
    # format
    if file_format is None:
        if input_file.endswith('.csv'):
            file_format = 'csv'
        elif input_file.endswith('.txt') or input_file.endswith('.tsv'):
            file_format = 'tab'
        elif input_file.endswith('.xls') or input_file.endswith('.xlsx'):
            file_format = 'excel'
    else:
        file_format = file_format.lower()
        
    # load via pandas IO
    if file_format == 'csv':
        df = pd.read_csv(input_file, sep=';', header=header)
    elif file_format == 'tab':
        df = pd.read_csv(input_file, sep='\t', header=header)
    elif file_format == 'excel':
        xls = pd.ExcelFile(input_file)
        df = pd.read_excel(xls, header=header)
    else:
        raise ValueError('Setup file not in usable format')
    
    # caps-invariant column IDs
    df.columns = [x.lower() for x in df.columns]
    
    # checking format of table
    check_df_setup(df)

    # filtering NaN for required columns
    df.dropna(subset=['sample type', 'sample labware name', 'sample location'],
              inplace=True)
    
    # making sure labware names are "TECAN worklist friendly"
    df = Utils.rm_special_chars(df, 'sample labware name')
    
    # assert & return
    assert df.shape[1] > 1, 'Input file is only 1 column; wrong delimiter used?'    
    return df


def check_df_setup(df_setup):
    """Assertions of df_conc object formatting
    """
    # checking for column IDs
    col_IDs = ('row', 'column', 'sample type', 
               'sample labware name', 'sample labware type',
               'sample location', 'sample volume',
               'mm name', 'mm volume', 'water volume')
    msg = 'Column "{}" not found (captilization invariant)'
    for x in col_IDs:
        if not x in df_setup.columns:
            raise ValueError(msg.format(x))

    # checking sample locations (>=1)
    msg = 'ERROR (SetupFile, line={}): location is < 1'
    for i,loc in enumerate(df_setup['sample location']):
        if loc < 1:
            print(msg.format(i), file=sys.stderr)
            
    # checking sample conc
    msg = 'ERROR (setupFile, line={}): volume is < 0'
    for i,vol in enumerate(df_setup['sample volume']):
        assert np.isnan(vol) or vol >= 0.0, msg.format(i)
    for i,vol in enumerate(df_setup['mm volume']):
        assert np.isnan(vol) or vol >= 0.0, msg.format(i)
    for i,vol in enumerate(df_setup['water volume']):
        assert np.isnan(vol) or vol >= 0.0, msg.format(i)

    # removing "tube" from end of labware type (if present)
    Utils.rm_tube(df_setup, 'sample labware type')


def check_rack_labels(df_setup):
    """Removing '.' for rack labels (causes execution failures)
    """
    cols = ['sample labware name', 'mm name', 'dest_labware_name']
    for x in cols:
        df_setup[x] = [str(y).replace('.', '_') for y in df_setup[x].tolist()]
    return df_setup
        
def plate2robot_loc(row_val, col_val, n_wells):
    """Changing positioning from row (letter) and column (number)
    to just numeric position (column-wise) on the plate,
    which is needed for TECAN robot positioning. 
    Using index for identifying well number in plate
    [args]
    row_val: string
    col_vol: string
    plate_type: string; plate type to determine well location indexing
    """    
    # index for converting row to numeric
    idx = string.ascii_uppercase
    idx = {x:i+1 for i,x in enumerate(idx)}    
    row_val = idx[row_val]

    # getting location on plate
    msg = 'Destination location "{}" is out of range'
    if n_wells == 96:
        loc = (col_val - 1) * 8 + row_val
        assert loc > 0 and loc <= 96, msg.format(loc)
    elif n_wells == 384:
        loc = (col_val - 1) * 16 + row_val
        assert loc > 0 and loc <= 384, msg.format(loc)
    else:
        msg = 'Number of wells is not valid: "{}"'
        raise ValueError(msg.format(plate_type))
    return loc


def add_dest(df_setup, dest_labware_name, dest_labware_type, n_wells=96):
    """Setting destination locations for samples & reagents
    Adding to df_conc:
      [dest_labware, dest_location]
    """
    # setting destination labware
    df_setup['dest_labware_name'] = dest_labware_name
    df_setup['dest_labware_type'] = dest_labware_type
    
    # setting destination location based on plate layout 
    func = lambda x: plate2robot_loc(x['row'], x['column'], n_wells=n_wells)
    df_setup['dest_target_position'] = df_setup.apply(func, 1)

def reorder_384well(df, reorder_col):
    """Reorder values so that the odd, then the even locations are
    transferred. This is faster for a 384-well plate
    df: pandas.DataFrame
    reorder_col: column name to reorder
    """
    df['TECAN_sort_IS_EVEN'] = [x % 2 == 0 for x in df[reorder_col]]
    df.sort_values(by=['TECAN_sort_IS_EVEN', reorder_col], inplace=True)
    df = df.drop('TECAN_sort_IS_EVEN', 1)
    df.index = range(df.shape[0])
    return df

def pip_mastermixes(df_setup, gwl, src_labware_type,
                    liq_cls='Mastermix Free Single',
                    n_tip_reuse=1):
    """Writing worklist commands for aliquoting mastermix.
    Re-using tips
    """
    # split by mastermix names (different master mixes)
    MM_names = np.unique(df_setup['mm name'])
    
    gwl.add(Fluent.Comment('Mastermixes'))
    for i,MM_name in enumerate(MM_names):
        # partitioning to just focal mastermix
        df = df_setup.loc[df_setup['mm name'] == MM_name]        
        df.reset_index(inplace=True)

        # all same volumes for mastermix?
        pip_mastermix(df, gwl=gwl,
                      MM_name=MM_name,                                     
                      src_labware_type=src_labware_type,
                      liq_cls=liq_cls,
                      n_tip_reuse=n_tip_reuse)


def pip_mastermix(df_map, gwl, MM_name, src_labware_type,
                  liq_cls='Mastermix Free Single',
                  n_tip_reuse=1):  
    """Dispense of particular mastermix
    """
    # df copy
    df = df_map.copy()
    ## ordering df for proper tip reuse
    x = cycle(range(8))
    df['CHANNEL_ORDER'] = [next(x) for y in range(df.shape[0])]
    x = cycle(range(n_tip_reuse))
    df['TIP_BATCH'] = Utils.tip_batch(df['CHANNEL_ORDER'], n_tip_reuse)
    df.sort_values(by=['TIP_BATCH',
                       'CHANNEL_ORDER',
                       'dest_target_position'], inplace=True)
    df.reset_index(inplace=True)
    
    # iterating mastermix records in setup table (single mastermix)
    gwl.add(Fluent.Comment('Mastermix: {}'.format(MM_name)))
    for i in range(df.shape[0]):
        # aspiration
        asp = Fluent.Aspirate()
        asp.RackLabel = '{0} MM[{1:0>3}]'.format(MM_name, 1)
        asp.RackType = src_labware_type
        asp.Position = 1
        asp.Volume = df.loc[i,'mm volume']
        asp.LiquidClass = liq_cls
        gwl.add(asp)
        
        # dispensing
        disp = Fluent.Dispense()
        disp.RackLabel = df.loc[i,'dest_labware_name']
        disp.RackType = df.loc[i,'dest_labware_type']
        disp.Position = df.loc[i,'dest_target_position']
        disp.Volume = df.loc[i,'mm volume']
        asp.LiquidClass = liq_cls
        gwl.add(disp)

        # waste
        if (i + 1) % n_tip_reuse == 0 or i + 1 == df.shape[0]:
            gwl.add(Fluent.Waste())
                
    # finish section
    gwl.add(Fluent.Break())

def pip_mastermix_multi_disp(df, gwl, MM_name, src_labware_type, multi_disp=6,
                             liq_cls='Mastermix Free Multi'):
    """Writing worklist commands for aliquoting mastermix.
    Re-using tips
    """
    # assertions
    cols = ['dest_labware_name', 'dest_labware_type']
    assert df.drop_duplicates(cols).shape[0] == 1
    
    # getting wells to exclude
    lw_type = df.loc[0,'dest_labware_type']
    n_wells = gwl.db.get_labware_wells(lw_type)
    all_wells = [x+1 for x in range(n_wells)]
    target_pos = df['dest_target_position'].tolist()
    to_exclude = set(all_wells) - set(target_pos)

    # creating reagnet distribution command
    rd = Fluent.Reagent_distribution()
    rd.SrcRackLabel = '{0} MM[{1:0>3}]'.format(MM_name, 1)
    rd.SrcRackType = '1.5ml Eppendorf waste'
    rd.SrcPosStart = 1
    rd.SrcPosEnd = 1
    # dispense parameters
    rd.DestRackLabel = df.loc[0,'dest_labware_name']
    rd.DestRackType = df.loc[0,'dest_labware_type']
    rd.DestPosStart = 1
    rd.DestPosEnd = n_wells
    # other
    rd.Volume = df.loc[0,'mm volume']
    rd.LiquidClass = liq_cls
    rd.NoOfDiTiReuses = 2
    rd.NoOfMultiDisp = multi_disp
    rd.Direction = 0
    rd.ExcludedDestWell = ';'.join([str(x) for x in list(to_exclude)])
    # adding to gwl object
    gwl.add(rd)
    
    # adding break
    gwl.add(Fluent.Break())        
        
def pip_samples(df_setup, gwl, liq_cls='Water Contact Wet Single'):
    """Commands for aliquoting samples into distination plate
    """

    gwl.add(Fluent.Comment('Samples'))    
    # filtering 'nan' from table
    x = pd.notnull(df_setup['sample labware name'])
    y = pd.notnull(df_setup['sample labware type'])
    z = pd.notnull(df_setup['sample location'])
    df = df_setup.loc[x & y & z]
    df.reset_index(inplace=True)
    if df.shape[0] < df_setup.shape[0]:
        msg = 'WARNING: some samples skipped due to missing values!'
        print(msg, file=sys.stderr)
    
    # for each Sample, create asp/dispense commands
    for i in range(df.shape[0]):
        # aspiration
        asp = Fluent.Aspirate()
        
        asp.RackLabel = df.loc[i,'sample labware name']
        asp.RackType = df.loc[i,'sample labware type']
        asp.Position = df.loc[i, 'sample location']
        asp.Volume = df.loc[i,'sample volume']
        asp.LiquidClass = liq_cls
        gwl.add(asp)
        
        # dispensing
        disp = Fluent.Dispense()
        disp.RackLabel = df.loc[i,'dest_labware_name']
        disp.RackType = df.loc[i,'dest_labware_type']
        disp.Position = df.loc[i,'dest_target_position']
        disp.Volume = df.loc[i,'sample volume']
        asp.LiquidClass = liq_cls
        gwl.add(disp)
        
        # waste (no tip re-use)
        gwl.add(Fluent.Waste())
        
    gwl.add(Fluent.Break())

        
def pip_water(df_setup, gwl, src_labware_type, 
              liq_cls='Water Contact Wet Single'):
    """Writing worklist commands for aliquoting water
    Using single asp-disp.
    """
    gwl.add(Fluent.Comment('Water')) 
    
    # filtering 'nan' from table
    x = pd.notnull(df_setup['water volume'])
    df = df_setup.loc[x]
    df.index = range(df.shape[0])
    if df.shape[0] < df_setup.shape[0]:
        msg = 'WARNING: water asp/disp for some samples skipped due to missing "water volume" values!'
        print(msg, file=sys.stderr)
    
    # for each Sample, create asp/dispense commands
    for i in range(df.shape[0]):
        if df.loc[i,'water volume'] <= 0:
            msg = 'WARNING: skipping water asp/disp for sample (volue <= 0)'
            print(msg, file=sys.stderr)
            continue
        
        # aspiration
        asp = Fluent.Aspirate()
        asp.RackLabel = 'Water source[{0:0>3}]'.format(1)        
        asp.RackType = src_labware_type
        asp.Position = 1
        asp.Volume = df.loc[i,'water volume']
        asp.LiquidClass = liq_cls
        gwl.add(asp)
        
        # dispensing
        disp = Fluent.Dispense()
        disp.RackLabel = df.loc[i,'dest_labware_name']
        disp.RackType = df.loc[i,'dest_labware_type']
        disp.Position = df.loc[i,'dest_target_position']
        disp.Volume = df.loc[i,'water volume']
        asp.LiquidClass = liq_cls
        gwl.add(disp)
        
        # waste (no tip re-use)
        gwl.add(Fluent.Waste())
        
    gwl.add(Fluent.Break())

def add_error(x, error_perc):
    if x is None:
        return None
    return x * (1.0 + error_perc / 100.0)

def write_report_line(outFH, subject, volume, round_digits=1, error_perc=None):
    if volume is None:
        v = 'NA'
    else:
        if error_perc is not None:
            volume = add_error(volume, error_perc)
        v = round(volume, round_digits)
    outFH.write('{}:\t{}\n'.format(subject, v))

def write_report(df_setup, MM_name, outFH):
    """Writing a report on the qPCR setup
    """
    # calculating total volumes
    n_rxn = df_setup.shape[0]
    ## total mastermix
    total_mm_volume = np.sum(df_setup['mm volume'])
    ## total water
    total_water_volume = np.sum(df_setup['water volume'])

    # report
    # number of samples
    outFH.write('# PCR REPORT\n')
    outFH.write('MasterMix: {}\n'.format(MM_name))
    outFH.write('Number of total rxns:\t{}\n'.format(n_rxn))
    ## raw total volumes
    outFH.write('# Total reagent volumes (ul)\n')
    write_report_line(outFH, 'MasterMix', total_mm_volume)
    write_report_line(outFH, 'Water', total_water_volume)
    ## with pipetting error
    outFH.write('# Total reagent volumes + 10% extra (ul)\n')
    write_report_line(outFH, 'MasterMix', total_mm_volume, error_perc=10)
    write_report_line(outFH, 'Water', total_water_volume, error_perc=10)
    ## end
    outFH.write('\n')

    
# main
if __name__ == '__main__':
    pass


