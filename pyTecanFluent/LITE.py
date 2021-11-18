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


# functions
def get_desc():
    desc = 'Create TECAN Fluent worklist instructions for LITE NGS library prep'
    return desc

def parse_args(test_args=None, subparsers=None):
    desc = get_desc()
    epi = """DESCRIPTION:
    Convert an input table of samples into a GWL file, which is used by the TECAN
    robot to conduct the NGS LITE (hacked Nextera) library prep.
    The extra columns in the mapping file designate the SOURCE of samples and primers;
    the DESTINATION (plate & well) is set by this script. 

    EXTRA COLUMNS in MAPPING FILE:

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
        parser = subparsers.add_parser('LITE', description=desc, epilog=epi,
                                       formatter_class=argparse.RawTextHelpFormatter)
    else:
        parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                         formatter_class=argparse.RawTextHelpFormatter)

    # args
    ## I/O
    groupIO = parser.add_argument_group('I/O')
    groupIO.add_argument('mapfile', metavar='MapFile', type=str,
                         help='A QIIME-formatted mapping file with extra columns (see below)')
    groupIO.add_argument('--rows', type=str, default='all',
                         help='Which rows of the mapping file to use (eg., "all"=all rows; "1-48"=rows1-48; "1,3,5-6"=rows1+3+5+6), (default: %(default)s)')
    groupIO.add_argument('--prefix', type=str, default='TECAN_LITE',
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
    rgnt = parser.add_argument_group('Reagents')
    rgnt.add_argument('--tag-mm-volume', type=float, default=3.0,
                      help='Tagmentation MasterMix volume per well (default: %(default)s)')
    rgnt.add_argument('--pcr-mm-volume', type=float, default=18.0,
                      help='PCR MasterMix volume per well (default: %(default)s)')
    rgnt.add_argument('--sample-volume', type=float, default=2.0,
                      help='Sample volume per PCR (default: %(default)s)')
    rgnt.add_argument('--primer-volume', type=float, default=2.0,
                         help='Primer volume per PCR, assuming foward+reverse are already combined (default: %(default)s)')
    rgnt.add_argument('--error-perc', type=float, default=10.0,
                        help='Percent of extra total reagent volume to include (default: %(default)s)')
    rgnt.add_argument('--tag-mm-labware-type', type=str, default='1.5ml Eppendorf waste',
                      help='Tagmentation: labware type for mastermix (default: %(default)s)')
    rgnt.add_argument('--pcr-mm-labware-type', type=str, default='25ml_1 waste',
                      help='PCR: labware type for mastermix (default: %(default)s)')
    
    # Liquid classes
    liq = parser.add_argument_group('Liquid classes')
    liq.add_argument('--tag-mm-liq', type=str, default='MasterMix Free Single Wall Disp',
                      help='Tagmentation: Mastermix liquid class (default: %(default)s)')
    liq.add_argument('--pcr-mm-liq', type=str, default='MasterMix Free Single',
                      help='PCR: Mastermix liquid class (default: %(default)s)')
    liq.add_argument('--sample-liq', type=str, default='Water Free Single Wall Disp',
                      help='Sample liquid class (default: %(default)s)')
    liq.add_argument('--primer-liq', type=str, default='Water Free Single Wall Disp',
                     help='Primer liquid class (default: %(default)s)')
    liq.add_argument('--tag-n-tip-reuse', type=int, default=4,
                     help='Tagmentation: number of tip reuses for multi-dispense (default: %(default)s)')
    liq.add_argument('--pcr-n-tip-reuse', type=int, default=4,
                     help='PCR: number of tip reuses for multi-dispense (default: %(default)s)')
    
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
    gwl_file, lw_file, report_file, df_file = main_PCR(df_map, args)
    
    # Return
    return gwl_file, lw_file, report_file, df_file

def main_tagmentation(df_map, args):
    """Tagmentation step of the LITE method
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

    # dispensing reagents (greater volume first)
    if args.tag_mm_volume <= 0 and args.sample_volume <= 0:
        return df_map
    elif args.tag_mm_volume >= args.sample_volume:
        ## mastermix
        pip_mastermix(df_map, gwl,
                      mm_labware_type=args.tag_mm_labware_type,
                      mm_volume=args.tag_mm_volume, 
                      liq_cls=args.tag_mm_liq,
                      n_tip_reuse=args.tag_n_tip_reuse)
        ## samples
        pip_samples(df_map, gwl, sample_volume=args.sample_volume, liq_cls=args.sample_liq)
    else:
        ## samples
        pip_samples(df_map, gwl, sample_volume=args.sample_volume, liq_cls=args.sample_liq)
        ## mastermix
        pip_mastermix(df_map, gwl,
                      mm_labware_type=args.tag_mm_labware_type,
                      mm_volume=args.tag_mm_volume, 
                      liq_cls=args.tag_mm_liq,
                      n_tip_reuse=args.tag_n_tip_reuse)

    ## writing out worklist (gwl) file
    gwl_file = args.prefix + '_tag.gwl'
    gwl.write(gwl_file)
    
    # Report (total volumes; sample truncation; samples)
    report_file = args.prefix + '_tag_report.txt'
    with open(report_file, 'w') as repFH:
        write_report(df_map, outFH=repFH,
                     mm_volume=args.tag_mm_volume,
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
    """PCR step of the LITE method
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
    pip_mastermix(df_map, gwl,
                  mm_labware_type=args.pcr_mm_labware_type,
                  mm_volume=args.pcr_mm_volume, 
                  liq_cls=args.pcr_mm_liq,
                  n_tip_reuse=args.pcr_n_tip_reuse)

    ## primers into tagmentation dest plate
    if args.primer_volume > 0:
        pip_primers(df_map, gwl,
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

    return gwl_file, lw_file, report_file, df_file
       
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

    # distination name
    args.dest_name = Utils.rm_special_chars(args.dest_name)
    
    # volumes
    ## mastermix
    if args.tag_mm_volume < 0:
        args.tag_mm_volume = 0
    if args.pcr_mm_volume < 0:
        args.pcr_mm_volume = 0        
    ## primers
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
    req_cols = ['TECAN_sample_labware_name', 'TECAN_sample_labware_type',
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
    df_map = Utils.rm_special_chars(df_map, 'TECAN_sample_labware_name')
    # removing "tube" from end of labware type (if present)
    Utils.rm_tube(df_map, 'TECAN_primer_labware_type')

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
                
def pip_samples(df_map, gwl, sample_volume, liq_cls='Water Free Single'):
    """Commands for aliquoting samples to each PCR rxn
    """
    gwl.add(Fluent.Comment('Samples'))
    # for each Sample-PCR, write out asp/dispense commands
    for i in range(df_map.shape[0]):
        if sample_volume <= 0:
            continue
        
        # aspiration
        asp = Fluent.Aspirate()
        asp.RackLabel = df_map.loc[i,'TECAN_sample_labware_name']
        asp.RackType = df_map.loc[i,'TECAN_sample_labware_type']
        asp.Position = df_map.loc[i,'TECAN_sample_target_position']
        asp.Volume = sample_volume
        asp.LiquidClass = liq_cls
        gwl.add(asp)

        # dispensing
        disp = Fluent.Dispense()
        disp.RackLabel = df_map.loc[i,'TECAN_dest_labware_name']
        disp.RackType = df_map.loc[i,'TECAN_dest_labware_type']
        disp.Position = df_map.loc[i,'TECAN_dest_target_position']
        disp.Volume = sample_volume
        disp.LiquidClass = liq_cls
        gwl.add(disp)

        # waste
        gwl.add(Fluent.Waste())
        
    # adding break
    gwl.add(Fluent.Break())

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
    if isinstance(df_biorad.index, pd.MultiIndex):
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


