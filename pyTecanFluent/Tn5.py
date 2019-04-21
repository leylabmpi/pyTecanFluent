from __future__ import print_function
# import
## batteries
import os
import sys
import argparse
import functools
from itertools import product,cycle
## 3rd party
import numpy as np
import pandas as pd
## package
from pyTecanFluent import Utils
from pyTecanFluent import Fluent
from pyTecanFluent import Labware
from pyTecanFluent import Tn5_pip


# functions
def get_desc():
    desc = 'Create TECAN Fluent worklist instructions for Tn5 library prep'
    return desc

def parse_args(test_args=None, subparsers=None):
    desc = get_desc()
    epi = """DESCRIPTION:
    Convert an input table of samples into a GWL file, which is used by the TECAN
    robot to conduct the NGS Tn5 (self-produced enzyme) library prep.
    The extra columns in the mapping file designate the SOURCE of samples and primers;
    the DESTINATION (plate & well) is set by this script. 

    The amount of Tn5 and buffer is determined based on the input. 

    Required columns in the input samples table:

    * "TECAN_sample_conc" = The sample concentrations (numeric value; units=ng/ul) 
    * "TECAN_sample_labware_name" = The sample labware name on the robot worktable. Whatever name you want to use! 
    * "TECAN_sample_labware_type" = The type of labware containing samples (eg., '96 Well Eppendorf TwinTec PCR')
    * "TECAN_sample_target_position" = The well or tube location (a number)
    * "TECAN_primer_labware_name" = The primer plate labware name on the robot worktable (eg., "515F-806R")
    * "TECAN_primer_labware_type" = The primer plate labware type on the robot worktable (eg., "96 Well Eppendorf TwinTec PCR"
    * "TECAN_primer_target_position" = The position (well) of your samples in your labware

    PRIMERS:
    * For single-indexed primers (eg., EMP primers), the non-barcoded primer should be pre-added to either
      the mastermix (adjust the volume used for this script!) or each of the barcoded primers.
    * If --primer-volume set to 0, then primers are skipped. 

    CONTROLS:
    * For the positive & negative controls, include them in the mapping file.
    * If the controls (or samples) are provided in a tube, use "1.5ml Eppendorf" 
      for the "TECAN_sample_labware_type" column. The labware name can be whatever you want.

    LABWARE:
    * In order to use a plate adapter, use "PCR Adapter 96 Well and 96 Well Eppendorf TwinTec PCR"

    MISC NOTES:
    * All volumes are in ul
    * Plate well locations are 1 to n-wells; numbering by column
    * PicoGreen should be added to the MasterMix *prior* to loading on robot
    """
    if subparsers:
        parser = subparsers.add_parser('Tn5', description=desc, epilog=epi,
                                       formatter_class=argparse.RawTextHelpFormatter)
    else:
        parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                         formatter_class=argparse.RawTextHelpFormatter)

    # args
    ## I/O
    groupIO = parser.add_argument_group('I/O')
    groupIO.add_argument('mapfile', metavar='MapFile', type=str,
                         help='A samples file with the required columns (see the help docs)')
    groupIO.add_argument('--rows', type=str, default='all',
                         help='Which rows of the sample file to use (eg., "all"=all rows; "1-48"=rows1-48; "1,3,5-6"=rows1+3+5+6), (default: %(default)s)')
    groupIO.add_argument('--prefix', type=str, default='TECAN_Tn5',
                         help='Output file name prefix (default: %(default)s)')

    ## Destination plate
    dest = parser.add_argument_group('Destination plate')
    dest.add_argument('--dest-name', type=str, default='Destination plate',
                      help='Distination labware name (default: %(default)s)')
    dest.add_argument('--dest-type', type=str, default='384 Well Biorad PCR',
                      choices=['96 Well Eppendorf TwinTec PCR',
                               'PCR Adapter 96 Well and 96 Well Eppendorf TwinTec PCR',
                               '384 Well Biorad PCR',
                               'PCR Adapter 384 Well and 384 Well Biorad PCR'],
                      help='Destination labware type (default: %(default)s)')
    dest.add_argument('--dest-start', type=int, default=1,
                      help='Start well number on destination plate (default: %(default)s)')
    
    ## Reagents
    tag_rgnt = parser.add_argument_group('Tagmentation Reagents')
    tag_rgnt.add_argument('--tag-rxn-volume', type=float, default=20.0,
                      help='Total tagmentation rxn volume per well (default: %(default)s)')
    tag_rgnt.add_argument('--min-sample-volume', type=float, default=1.0,
                      help='Min sample volume per tag. rxn (default: %(default)s)')
    tag_rgnt.add_argument('--max-sample-volume', type=float, default=10.0,
                      help='Max sample volume per tag. rxn (default: %(default)s)')
    tag_rgnt.add_argument('--target-ng-sample', type=float, default=5.0,
                      help='Target amount of DNA (ng) to use per sample (default: %(default)s)')
    tag_rgnt.add_argument('--tag-Tn5-labware-type', type=str, default='2ml Eppendorf waste',
                      help='Labware type for Tn5 enzyme (default: %(default)s)')
    tag_rgnt.add_argument('--tag-buffer-labware-type', type=str, default='2ml Eppendorf waste',
                      help='Labware type for Tn5 buffer (default: %(default)s)')
    tag_rgnt.add_argument('--tag-n-tip-reuse', type=int, default=4,
                      help='Number of tip reuses for multi-dispense (default: %(default)s)')


    pcr_rgnt = parser.add_argument_group('PCR Reagents')    
    pcr_rgnt.add_argument('--pcr-mm-volume', type=float, default=24.0,
                          help='PCR MasterMix volume per well (default: %(default)s)')
    pcr_rgnt.add_argument('--primer-volume', type=float, default=6.0,
                          help='Primer volume per PCR, assuming foward+reverse are already combined (default: %(default)s)')
    pcr_rgnt.add_argument('--error-perc', type=float, default=10.0,
                          help='Percent of extra total reagent volume to include (default: %(default)s)')
    pcr_rgnt.add_argument('--pcr-mm-labware-type', type=str, default='25ml_1 waste',
                          help='Labware type for mastermix (default: %(default)s)')
    pcr_rgnt.add_argument('--pcr-n-tip-reuse', type=int, default=4,
                          help='PCR: number of tip reuses for multi-dispense (default: %(default)s)')
    
    # Liquid classes
    liq = parser.add_argument_group('Liquid classes')
    liq.add_argument('--tag-Tn5-liq', type=str, default='MasterMix Free Single Wall Disp',
                      help='Tagmentation: Mastermix liquid class (default: %(default)s)')
    liq.add_argument('--tag-buffer-liq', type=str, default='MasterMix Free Single Wall Disp',
                      help='Tagmentation: buffer liquid class (default: %(default)s)')
    liq.add_argument('--sample-liq', type=str, default='Water Free Single Wall Disp',
                      help='Sample liquid class (default: %(default)s)')    
    liq.add_argument('--pcr-mm-liq', type=str, default='MasterMix Free Single',
                      help='PCR: Mastermix liquid class (default: %(default)s)')
    liq.add_argument('--primer-liq', type=str, default='Water Free Single Wall Disp',
                     help='Primer liquid class (default: %(default)s)')
    liq.add_argument('--water-liq', type=str, default='Water Free Single Wall Disp',
                      help='Water liquid class (default: %(default)s)')    

    misc = parser.add_argument_group('Misc')     
    misc.add_argument('--water-labware-type', type=str, default='25ml_1 waste',
                      help='Labware type for water (default: %(default)s)')
    
    # running test args
    if test_args:
        args = parser.parse_args(test_args)
        return args

    return parser

def main(args=None):
    # Input
    if args is None:
        args = parse_args()
    check_args(args)
    
    # Import
    df_map = map2df(args.mapfile, args, row_select=args.rows)
    
    # Making destination dataframe
    df_map = add_dest(df_map,
                      dest_labware=args.dest_name,
                      dest_type=args.dest_type,
                      dest_start=args.dest_start)
    df_map = check_rack_labels(df_map)
    
    # tagmentation assay setup
    df_map = main_tagmentation(df_map, args)
    
    # PCR assay setup
    main_PCR(df_map, args)
    
    # Return
    return None

def calc_Tn5_volume(x, min_vol = 0.01):
    """Calculating the volume in ul of Tn5 (y)
    base on DNA conc. (ng) input (x)
    """
    y = None
    if x >= 0 and x < 3:
        y = x * 0.1678 + 0.05
    elif x >= 3 and x < 15:
        y = x * 0.02 + 0.49
    elif x >= 15 and x < 25:
        y = x * 0.227 - 2.61
    elif x >= 25:
        y = x * 0.12
    else:
        raise ValueError('Logic error')
    
    if y < min_vol:
        y = min_vol
        
    return y 

def calc_Tn5_buffer_volume(x, min_vol = 1.0):
    """Calculating the amount of buffer to use depending on
    the Tn5 volume
    """
    y = x * 4.0 / 3.0
    if y < min_vol:
        y = min_vol
    return y

def calc_Tn5_H2O_volume(DNA_vol, Tn5_vol, buf_vol, total_vol):
    y = total_vol - (DNA_vol + Tn5_vol + buf_vol)
    if y < 0:
        raise ValueError('H2O volume is < 0')
    return y

def set_Tn5_volume(Tn5_vol):
    """calculating the volume for the necessary 10-fold dilution 
    """
    if Tn5_vol >= 1.0:
        return Tn5_vol
    elif Tn5_vol < 1.0 and Tn5_vol >= 0.1:
        # 1:10 dilution
        return Tn5_vol * 10
    elif Tn5_vol < 0.1 and Tn5_vol >= 0.01:
        # 1:100 dilution
        return Tn5_vol * 100
    else:
        msg = 'Calculated undiluted Tn5 volume too low: {}'
        raise ValueError(msg.format(Tn5_vol))

def set_Tn5_dilution_volume(DNA_ng, dilution=1):
    """calculating the volume for the necessary 10-fold dilution 
    """
    # undiluted Tn5 volume
    Tn5_vol = calc_Tn5_volume(DNA_ng)

    # volume needed if diluting
    if dilution == 1:
        if Tn5_vol >= 1.0:
            return Tn5_vol
        else:
            return np.nan
    elif dilution == 10:
        if Tn5_vol < 1.0 and Tn5_vol >= 0.1:
            return Tn5_vol * 10
        else:
            return np.nan
    elif dilution == 100:
        if Tn5_vol < 0.1 and Tn5_vol >= 0.01:
            return Tn5_vol * 100
        else:
            return np.nan
    else:
        msg = 'Calculated undiluted Tn5 volume too low: {}'
        raise ValueError(msg.format(Tn5_vol))
    
def calc_DNA_volume(DNA_conc, target_DNA_ng=5.0,
                    min_DNA_vol=1.0, max_DNA_vol=10.0,
                    total_rxn_vol=20.0):
    """DNA conc. in ng/ul
    """
    # setting rxn_volume-limited max_DNA_vol
    ul_step = -0.5
    attempt_ul_DNA = max_DNA_vol
    while 1:
        if attempt_ul_DNA <= min_DNA_vol:
            # lowest possible volume reached
            max_DNA_vol = min_DNA_vol
            break
        
        # calculating reagent volumes
        attempt_ng_DNA = DNA_conc * attempt_ul_DNA
        attempt_ul_Tn5 = calc_Tn5_volume(attempt_ng_DNA)
        attempt_ul_buf = calc_Tn5_buffer_volume(attempt_ul_Tn5)
        attempt_ul_Tn5 = set_Tn5_volume(attempt_ul_Tn5)
        attempt_ul_H2O = total_rxn_vol - (attempt_ul_DNA + attempt_ul_Tn5 + attempt_ul_buf)
        
        # checking values
        if attempt_ul_H2O < 0:
            if attempt_ul_DNA == min_DNA_vol:
                msg = 'Min sample volume set too high for total rxn volume!'
                raise ValueError(msg)
            # too much (DNA + reagents) for total rxn volume
            attempt_ul_DNA += ul_step
        else:
            # found 'true' max volume
            max_DNA_vol = attempt_ul_DNA
            break

    #print('Max DNA volume: {}'.format(max_DNA_vol))

    # determining amount of DNA volume to use
    ul_step = 0.5
    max_vol_exceeded = False
    attempt_ul_DNA = min_DNA_vol
    # iterating; trying to get to target sample conc
    while 1:
        attempt_ng_DNA = DNA_conc * attempt_ul_DNA
        attempt_ul_Tn5 = calc_Tn5_volume(attempt_ng_DNA)
        attempt_ul_buf = calc_Tn5_buffer_volume(attempt_ul_Tn5)
        attempt_ul_Tn5 = set_Tn5_volume(attempt_ul_Tn5)
        attempt_ul_H2O = total_rxn_vol - (attempt_ul_DNA + attempt_ul_Tn5 + attempt_ul_buf)

        # hitting constraints
        if attempt_ul_DNA > max_DNA_vol or attempt_ul_H2O < 0:
            # too much reagents volume
            ## iterate backwards and select first that fits constraint
            if attempt_ul_DNA == min_DNA_vol:
                raise ValueError('logic error')
            attempt_ul_DNA -= ul_step
            if attempt_ul_DNA < min_DNA_vol:
                return min_DNA_vol
            max_vol_exceeded = True
            continue
        elif max_vol_exceeded is True:
            return attempt_ul_DNA    
        
        # checking 
        if attempt_ng_DNA < target_DNA_ng:
            # try to increase volume
            attempt_ul_DNA += ul_step
        elif attempt_ng_DNA >= target_DNA_ng:
            # value reached
            return attempt_ul_DNA
        else:
            raise ValueError('logic error')    

def calc_tag_volumes(df_map, args):
    """Determining the volume to Tn5, buffer, and water based
    on DNA input and total MasterMix volume
    """
    # calculating volume of DNA to use
    func = functools.partial(calc_DNA_volume,
                             target_DNA_ng = args.target_ng_sample,
                             min_DNA_vol = args.min_sample_volume,
                             max_DNA_vol = args.max_sample_volume,
                             total_rxn_vol = args.tag_rxn_volume)
    df_map['TECAN_sample_ul'] = df_map['TECAN_sample_conc'].apply(func)

    cols = ['TECAN_sample_ul', 'TECAN_sample_conc']
    df_map['TECAN_sample_ng'] = df_map[cols].apply(lambda x: x[0] * x[1], axis=1)
    
    # calculating Tn5
    df_map['TECAN_Tn5_ul'] = df_map['TECAN_sample_ng'].apply(calc_Tn5_volume)
    # calculating Tn5 buffer
    df_map['TECAN_Tn5_buffer_ul'] = df_map['TECAN_Tn5_ul'].apply(calc_Tn5_buffer_volume)
    # calculating true Tn5
    df_map['TECAN_Tn5_ul'] = df_map['TECAN_Tn5_ul'].apply(set_Tn5_volume)

    # calculating Water
    func = functools.partial(calc_Tn5_H2O_volume, total_vol = args.tag_rxn_volume)
    cols = ['TECAN_sample_ul', 'TECAN_Tn5_ul', 'TECAN_Tn5_buffer_ul']
    df_map['TECAN_Tn5_H2O_ul'] = df_map[cols].apply(lambda x: func(x[0], x[1], x[2]), axis=1)

    # Tn5 dilutions
    ## no dilution
    func = functools.partial(set_Tn5_dilution_volume, dilution=1)
    df_map['TECAN_Tn5_1fd_ul'] = df_map['TECAN_sample_ng'].apply(func)
    ## 1:10 dilution
    func = functools.partial(set_Tn5_dilution_volume, dilution=10)
    df_map['TECAN_Tn5_10fd_ul'] = df_map['TECAN_sample_ng'].apply(func)
    ## 1:100 dilution
    func = functools.partial(set_Tn5_dilution_volume, dilution=100)
    df_map['TECAN_Tn5_100fd_ul'] = df_map['TECAN_sample_ng'].apply(func)

    # rounding values
    idx = {'TECAN_Tn5_ul' : 2, 'TECAN_Tn5_ul' : 2,
           'TECAN_Tn5_buffer' : 2, 'TECAN_Tn5_H2O_ul' : 2,
           'TECAN_Tn5_1fd_ul' : 2, 'TECAN_Tn5_10fd_ul' : 2,
           'TECAN_Tn5_100fd_ul' : 2}
    df_map = df_map.round(idx)
    
    # ret
    return df_map
    
def main_tagmentation(df_map, args):
    """Tagmentation step of the Tn5 method
    """
    # calculating volumes
    df_map = calc_tag_volumes(df_map, args)
    
    # gwl construction
    TipTypes = ['FCA, 1000ul SBS', 'FCA, 200ul SBS',
                'FCA, 50ul SBS', 'FCA, 10ul SBS']     
    gwl = Fluent.gwl(TipTypes)

    # Reordering dest if plate type is 384-well
    df_map = Utils.reorder_384well(df_map, gwl,
                                   labware_name_col='TECAN_dest_labware_name',
                                   labware_type_col='TECAN_dest_labware_type',
                                   position_col='TECAN_dest_target_position')


    
    # dispensing reagents (greater volume first)
    ## Tn5 buffer
    Tn5_pip.pip_Tn5_buffer(df_map, gwl,
                           src_labware_type=args.tag_buffer_labware_type,
                           liq_cls=args.tag_buffer_liq,
                           n_tip_reuse=args.tag_n_tip_reuse)
    ## Tn5
    Tn5_pip.pip_Tn5(df_map, gwl,
                    src_labware_type=args.tag_Tn5_labware_type,
                    liq_cls=args.tag_Tn5_liq,
                    n_tip_reuse=args.tag_n_tip_reuse)
    ## water
    Tn5_pip.pip_tag_water(df_map, gwl,
                          src_labware_type=args.water_labware_type,
                          liq_cls=args.water_liq,
                          n_tip_reuse=args.tag_n_tip_reuse)
    ## samples
    Tn5_pip.pip_samples(df_map, gwl,
                        liq_cls=args.sample_liq)
    
    ## writing out worklist (gwl) file
    gwl_file = args.prefix + '_tag.gwl'
    gwl.write(gwl_file)
    
    # Report (total volumes; sample truncation; samples)
    report_file = args.prefix + '_tag_report.txt'
    with open(report_file, 'w') as repFH:
        write_report(df_map, outFH=repFH,
                     mm_volume=args.tag_rxn_volume,
                     error_perc=args.error_perc)
        
    # making labware table
    lw = Labware.labware()
    lw.add_gwl(gwl)
    lw_df = lw.table()
    lw_file = args.prefix + '_tag_labware.txt'
    lw_df.to_csv(lw_file, sep='\t', index=False)
        
    # Mapping file with destinations
    df_file = args.prefix + '_tag_map.txt'
    df_map['TECAN_dest_target_position'] = df_map['TECAN_dest_target_position'].astype(int)
    df_map.to_csv(df_file, sep='\t', index=False, na_rep='NA')
    
    # status on files written
    Utils.file_written(gwl_file)
    Utils.file_written(lw_file)
    Utils.file_written(report_file)
    Utils.file_written(df_file)

    # returning modified df_map
    return df_map

def main_PCR(df_map, args):
    """PCR step of the Tn5 method
    """
    # gwl construction
    TipTypes = ['FCA, 1000ul SBS', 'FCA, 200ul SBS',
                'FCA, 50ul SBS', 'FCA, 10ul SBS']     
    gwl = Fluent.gwl(TipTypes)
    
    # Reordering dest if plate type is 384-well
    df_map = Utils.reorder_384well(df_map, gwl,
                                   labware_name_col='TECAN_dest_labware_name',
                                   labware_type_col='TECAN_dest_labware_type',
                                   position_col='TECAN_dest_target_position')

    ## mastermix into tagmentation dest plate
    Tn5_pip.pip_mastermix(df_map, gwl,
                          mm_labware_type=args.pcr_mm_labware_type,
                          mm_volume=args.pcr_mm_volume, 
                          liq_cls=args.pcr_mm_liq,
                          n_tip_reuse=args.pcr_n_tip_reuse)

    ## primers into tagmentation dest plate
    if args.primer_volume > 0:
        Tn5_pip.pip_primers(df_map, gwl,
                            prm_volume=args.primer_volume,
                            liq_cls=args.primer_liq)
    
    ## writing out worklist (gwl) file
    gwl_file = args.prefix + '_pcr.gwl'
    gwl.write(gwl_file)
    
    # Report (total volumes; sample truncation; samples)
    report_file = args.prefix + '_pcr_report.txt'
    with open(report_file, 'w') as repFH:
        write_report(df_map, outFH=repFH,
                     mm_volume=args.pcr_mm_volume,
                     prm_volume=args.primer_volume,
                     error_perc=args.error_perc)
        
    # making labware table
    lw = Labware.labware()
    lw.add_gwl(gwl)
    lw_df = lw.table()
    lw_file = args.prefix + '_pcr_labware.txt'
    lw_df.to_csv(lw_file, sep='\t', index=False)
        
    # Mapping file with destinations
    df_file = args.prefix + '_pcr_map.txt'
    df_map['TECAN_dest_target_position'] = df_map['TECAN_dest_target_position'].astype(int)
    df_map.to_csv(df_file, sep='\t', index=False, na_rep='NA')

    # Plate map file for Bio-Rad PrimePCR software (designates: sampleID <--> wellID)
    biorad_files = PrimerPCR_plate_map(df_map, prefix=args.prefix)
    
    # status on files written
    Utils.file_written(gwl_file)
    Utils.file_written(lw_file)
    Utils.file_written(report_file)
    Utils.file_written(df_file)
    for F in biorad_files:
        Utils.file_written(F)        
       
def check_args(args):
    """Checking user input
    """
    # destination start
    db = Fluent.db()
    db.get_labware(args.dest_type)
    wells = db.get_labware_wells(args.dest_type)
    if args.dest_start < 1 or args.dest_start > wells:
        msg = 'Destination start well # must be in range: 1-{}'
        raise ValueError(msg.format(wells))

    # rows in mapping file
    args.rows = Utils.make_range(args.rows, set_zero_index=True)
    
    # volumes
    ## Tn5
    if args.min_sample_volume < 0:
        args.min_sample_volume = 0
    if args.target_ng_sample < 0:
        args.target_ng_sample = 0
    ## PCR
    ### mastermix
    if args.tag_rxn_volume < 0:
        args.tag_rxn_volume = 0
    if args.pcr_mm_volume < 0:
        args.pcr_mm_volume = 0        
    ### primers
    if args.primer_volume < 0:
        args.primer_volume = 0
    
        
def map2df(mapfile, args, row_select=None):
    """Loading a mapping file as a pandas dataframe
    mapfile: string; mapping file path
    row_select: select particular rows of table in map file
    """
    # load via pandas IO
    if mapfile.endswith('.txt') or mapfile.endswith('.tsv'):
        df = pd.read_csv(mapfile, sep='\t')
    elif mapfile.endswith('.csv'):
        df = pd.read_csv(mapfile, sep=',')
    elif mapfile.endswith('.xls') or mapfile.endswith('.xlsx'):
        xls = pd.ExcelFile(mapfile)
        df = pd.read_excel(xls)
    else:
        raise ValueError('Mapping file not in usable format')

    # selecting particular rows
    if row_select is not None:
        df = df.iloc[row_select]

    # checking format
    df = check_df_map(df, args)
        
    # return
    return df

def missing_cols(df, req_cols):
    msg = 'Required column "{}" not found'
    for req_col in req_cols:
        if req_col not in df.columns.values:
            raise ValueError(msg.format(req_col))    

def filter_nan_cols(df, col):
    """Filtering nan columns if present
    """
    if df[col].isnull().values.any():
        msg = 'WARNING: nan values found in column: "{}". Filtering these nan rows\n'
        sys.stderr.write(msg.format(col))
        
    return df.dropna(subset=[col])
    
def check_df_map(df_map, args):
    """Assertions of df_map object formatting
    * Assumes `sample` field = 1st column
    df_map: map dataframe
    """
    # checking columns
    ## universal
    to_rename = {}
    if 'TECAN_labware_name' in df_map.columns.values and 'TECAN_sample_labware_name' not in df_map.columns.values:
        to_rename['TECAN_labware_name'] = 'TECAN_sample_labware_name'
    if 'TECAN_labware_type' in df_map.columns.values and 'TECAN_sample_labware_type' not in df_map.columns.values:
        to_rename['TECAN_labware_type'] = 'TECAN_sample_labware_type'
    if 'TECAN_target_position' in df_map.columns.values and 'TECAN_sample_target_position' not in df_map.columns.values:
        to_rename['TECAN_target_position'] = 'TECAN_sample_target_position'
    df_map.rename(columns=to_rename, inplace=True)    
    req_cols = ['TECAN_sample_conc',
                'TECAN_sample_labware_name', 'TECAN_sample_labware_type',
                'TECAN_sample_target_position', 'TECAN_primer_labware_name',
                'TECAN_primer_labware_type', 'TECAN_primer_target_position']
    missing_cols(df_map, req_cols)

    # adding samples column if not present
    if '#SampleID' in df_map.columns.values:
        df_map.rename(columns={'#SampleID':'SampleID'}, inplace=True)
    if 'SampleID' not in df_map.columns.values:
        df_map['SampleID'] = [x + 1 for x in range(df_map.shape[0])]

    # filtering nan columns if present
    for x in req_cols:
        df_map = filter_nan_cols(df_map, x)
        
    # checking for unique samples
    if any(df_map.duplicated('SampleID')):
        msg = 'WARNING: Duplicated sample values in the mapping file. Adding plate name to SampleID to make values unique'
        print(msg, file=sys.stderr)
        df_map['SampleID'] = df_map[['SampleID', 'TECAN_sample_labware_name']].apply(lambda x: '__'.join([str(y) for y in x]), axis=1)
        if any(df_map.duplicated('SampleID')):
            msg = 'ERROR: Duplicated sample values in the mapping file on the same plate! Not acceptable!'
            raise IOError(msg)

    # checking for unique barcode locations (NaN's filtered out)
    dups = df_map[req_cols].dropna().duplicated(keep=False)
    if any(dups):
        msg = 'WARNING: Duplicated barcodes in the mapping file!'
        print(msg, file=sys.stderr)        

    # making sure labware names are "TECAN worklist friendly"
    Utils.df_rm_special_chars(df_map, 'TECAN_sample_labware_name')

    # return
    return df_map

        
def check_rack_labels(df_map):
    """Removing '.' for rack labels (causes execution failures)
    """
    cols = ['TECAN_sample_labware_name',
            'TECAN_primer_labware_name',
            'TECAN_dest_labware_name']
    for x in cols:
        df_map[x] = [y.replace('.', '_') for y in df_map[x].tolist()]
    return df_map
        
def add_dest(df_map, dest_labware, dest_type='384 Well Biorad PCR', dest_start=1):
    """Setting destination locations for samples & primers.
    Making a new dataframe with:
    [dest_labware, dest_location]
    Joining to df_map
    """
    db = Fluent.db()    
    dest_start= int(dest_start)
    
    # labware type found in DB?
    db.get_labware(dest_type)
    positions = db.get_labware_wells(dest_type)

    # init destination df
    sample_col = 'SampleID'
    cols = [sample_col, 
            'TECAN_dest_labware_name',
            'TECAN_dest_labware_type', 
            'TECAN_dest_target_position']    
    ncol = len(cols)
    nrow = df_map.shape[0] 
    df_dest = pd.DataFrame(np.nan, index=range(nrow), columns=cols)

    # number of destination plates required
    n_dest_plates = round(df_dest.shape[0] / positions + 0.5, 0)
    if n_dest_plates > 1:
        msg = ('WARNING: Not enough wells for the number of samples.' 
        ' Using multiple destination plates')
        print(msg, file=sys.stderr)
            
    # filling destination df
    orig_dest_labware = dest_labware
    for i,sample in enumerate(df_map['SampleID']):
        # dest location
        dest_position = i + dest_start
        dest_position = positions if dest_position % positions == 0 else dest_position % positions 

        # destination plate name
        if n_dest_plates > 1:
            x = round((i + dest_start) / positions + 0.499, 0)        
            dest_labware = '{} {}'.format(orig_dest_labware, int(x))

        # adding values DF
        df_dest.iloc[i] = [sample, dest_labware, dest_type, dest_position]

    # df join (map + destination)
    assert df_map.shape[0] == df_dest.shape[0], 'df_map and df_dest are different lengths' 
    df_j = pd.merge(df_map, df_dest, on=sample_col, how='inner')
    assert df_j.shape[0] == df_dest.shape[0], 'map-dest DF join error'
    assert df_j.shape[0] > 0, 'DF has len=0 after adding destinations'

    # return
    return df_j

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

def PrimerPCR_plate_map(df_map, prefix='PrimerPCR', sep=','):
    """Create a PrimerPCR plate map table.
    This table is imported by PrimerPCR to designate: sampleID <--> wellID
    Table columns: Row     Column  *Target Name    *Sample Name
      "Row" : plate row
      "Column" : plate column
      "*Target Name" : Not needed (NA)
      "*Sample Name" : Sample in well
    Return: list of pandas dataframes (1 dataframe per file)
    """
    # Plate import file for Bio-Rad PrimePCR software
    dest_pos_max = 96 if df_map['TECAN_dest_target_position'].max() <= 96 else 384
    f = functools.partial(map2biorad, positions=dest_pos_max)
    df_biorad = df_map.groupby('TECAN_dest_labware_name').apply(f)
    biorad_files = []
    if isinstance(df_biorad.index, pd.core.index.MultiIndex):
        for labware in df_biorad.index.get_level_values(0).unique():
            biorad_file = prefix + '_BIORAD-{}.txt'.format(labware.replace(' ', '_'))
            df_biorad.loc[labware].to_csv(biorad_file, sep=sep, index=False, na_rep='')
            biorad_files.append(biorad_file)
    else:
        biorad_file = prefix + '_BIORAD.txt'
        df_biorad.to_csv(biorad_file, sep=sep, index=False, na_rep='')
        biorad_files.append(biorad_file)
    
    return biorad_files
        
def map2biorad(df_map, positions):
    """Making a table for loading into the Bio-Rad PrimePCR software
    columns are:  Row,Column,*Target Name,*Sample Name
    positions = max number of positions on qPCR plate
    Notes:
      Row = letter format
      Column = numeric
    """
    df_biorad = df_map.copy()
    df_biorad = df_biorad[['SampleID', 'TECAN_dest_target_position']]
    df_biorad.columns = ['*Sample Name', 'TECAN_dest_target_position']
    lw_utils = Labware.utils()
    f = functools.partial(lw_utils.position2well, wells = positions)
    x = df_map['TECAN_dest_target_position'].apply(f)
    x = pd.DataFrame([y for y in x], columns=['Row', 'Column'])
    df_biorad['Row'] = x.iloc[:,0].tolist()
    df_biorad['Column'] = x.iloc[:,1].tolist()

    df_biorad['*Target Name'] = np.nan
    df_biorad = df_biorad[['Row', 'Column', '*Target Name', '*Sample Name']]
    
    return df_biorad
        
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

def write_report(df_map, outFH, mm_volume=0.0,
                 prm_volume=0.0, error_perc=10.0):
    """Writing a report on
    """
    # calculating total volumes
    n_rxn = df_map.shape[0]
    ## total mastermix
    total_mm_volume = mm_volume * n_rxn

    # report
    # number of samples
    outFH.write('# RXN REPORT\n')
    outFH.write('Number of total RXNs:\t{}\n'.format(n_rxn))
    ## rxn volumes
    outFH.write('# Volumes per RXN (ul)\n')
    write_report_line(outFH, 'Master Mix', mm_volume)
    write_report_line(outFH, 'Primers', prm_volume)
    ## raw total volumes
    outFH.write('# Total volumes (ul)\n')
    write_report_line(outFH, 'Master Mix', total_mm_volume)
    ## total volumes with error
    outFH.write('# Total volumes with (ul; {}% error)\n'.format(error_perc))
    write_report_line(outFH, 'Master Mix', total_mm_volume, error_perc=error_perc)
    # samples
    outFH.write('')


# main
if __name__ == '__main__':
    pass


