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
from pyTecanFluent import Tn5_onBead_pip


# functions
def get_desc():
    desc = 'Create TECAN Fluent worklist instructions for Tn5 library prep'
    return desc

def parse_args(test_args=None, subparsers=None):
    desc = get_desc()
    epi = """DESCRIPTION:
    Convert an input table of samples into a GWL file, which is used by the TECAN
    robot to conduct the NGS Tn5 (self-produced enzyme) library prep.
    The extra columns in the mapping file designate the SOURCE of samples and primers.
    The DESTINATION (plate & well) is set by this script. 

    Required columns in the input samples table:
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
        parser = subparsers.add_parser('Tn5_onBead', description=desc, epilog=epi,
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
    groupIO.add_argument('--prefix', type=str, default='TECAN_Tn5-on-Bead',
                         help='Output file name prefix (default: %(default)s)')
    
    ## Reagents
    pcr_rgnt = parser.add_argument_group('PCR Reagents')
    pcr_rgnt.add_argument('--sup-volume', type=float, default=120.0,
                          help='Supernatant volume (default: %(default)s)')
    pcr_rgnt.add_argument('--mm-volume', type=float, default=23.0,
                          help='PCR MasterMix volume per well (default: %(default)s)')
    pcr_rgnt.add_argument('--primer-volume', type=float, default=6.0,
                          help='Primer volume per PCR, assuming foward+reverse are already combined (default: %(default)s)')
    pcr_rgnt.add_argument('--mm-labware-type', type=str, default='25ml_1 waste',
                          help='Labware type for mastermix (default: %(default)s)')
    
    # Liquid classes
    liq = parser.add_argument_group('Liquid classes')
    liq.add_argument('--sup-liq', type=str, default='Tn5-on-bead Supernatant Free Single',
                      help='Supernatant removal liquid class (default: %(default)s)')    
    liq.add_argument('--mm-liq', type=str, default='MasterMix Free Single',
                      help='PCR: Mastermix liquid class (default: %(default)s)')
    liq.add_argument('--primer-liq', type=str, default='Water Contact Wet Single',
                     help='Primer liquid class (default: %(default)s)')

    misc = parser.add_argument_group('Misc')     
    misc.add_argument('--error-perc', type=float, default=10.0,
                      help='Percent of extra total reagent volume to include (default: %(default)s)')
    
    # running test args
    if test_args:
        args = parser.parse_args(test_args)
        return args

    return parser

def main(args=None):
    """
    Main iterface
    """
    # Input
    if args is None:
        args = parse_args()
    check_args(args)
    
    # Import
    df_map = map2df(args.mapfile, args, row_select=args.rows)
    
    # Making destination dataframe
    df_map = add_dest(df_map)
    df_map = check_rack_labels(df_map)
    
    # PCR assay setup
    df_map, pcr_files = main_PCR(df_map, args)
    
    # Return
    return pcr_files
    
def main_PCR(df_map, args):
    """
    PCR step of the Tn5 method
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

    ## PCR mastermix into tagmentation dest plate
    Tn5_onBead_pip.pip_mastermix(df_map, gwl,
                                 mm_labware_type=args.mm_labware_type,
                                 mm_volume=args.mm_volume, 
                                 liq_cls=args.mm_liq,
                                 sup_volume=args.sup_volume,
                                 sup_rm_liq_cls=args.sup_liq)

    ## primers into tagmentation dest plate
    if args.primer_volume > 0:
        Tn5_onBead_pip.pip_primers(df_map, gwl,
                                   prm_volume=args.primer_volume,
                                   liq_cls=args.primer_liq)
    
    ## writing out worklist (gwl) file
    gwl_file = args.prefix + '_pcr.gwl'
    gwl.write(gwl_file)
    
    # Report (total volumes; sample truncation; samples)
    report_file = args.prefix + '_pcr_report.txt'
    with open(report_file, 'w') as repFH:
        write_pcr_report(df_map, outFH=repFH,
                         mm_volume=args.mm_volume,
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

    return df_map, [gwl_file, lw_file, report_file, df_file]
       
def check_args(args):
    """
    Checking user input
    """
    # destination start
    db = Fluent.db()
    #db.get_labware(args.dest_type)
    #wells = db.get_labware_wells(args.dest_type)
    #if args.dest_start < 1 or args.dest_start > wells:
    #    msg = 'Destination start well # must be in range: 1-{}'
    #    raise ValueError(msg.format(wells))

    # rows in mapping file
    args.rows = Utils.make_range(args.rows, set_zero_index=True)
    
    # volumes
    ## PCR
    ### mastermix
    if args.mm_volume < 0:
        args.mm_volume = 0        
    ### primers
    if args.primer_volume < 0:
        args.primer_volume = 0    
        
def map2df(mapfile, args, row_select=None):
    """
    Loading a mapping file as a pandas dataframe
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
    """
    Filtering nan columns if present
    """
    if df[col].isnull().values.any():
        msg = 'WARNING: nan values found in column: "{}". Filtering these nan rows\n'
        sys.stderr.write(msg.format(col))
        
    return df.dropna(subset=[col])
    
def check_df_map(df_map, args):
    """
    Assertions of df_map object formatting
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
    df_map = Utils.rm_special_chars(df_map, 'TECAN_primer_labware_name')

    # removing "tube" from end of labware type (if present)
    Utils.rm_tube(df_map, 'TECAN_sample_labware_type')
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
        
def add_dest(df_map):#, dest_labware, dest_type='384 Well Biorad PCR', dest_start=1):
    """
    Setting destination locations for samples & primers.
    Making a new dataframe with:
    [dest_labware, dest_location]
    Joining to df_map
    """
    df_map['TECAN_dest_labware_name'] = df_map['TECAN_sample_labware_name']
    df_map['TECAN_dest_labware_type'] = df_map['TECAN_sample_labware_type']
    df_map['TECAN_dest_target_position'] = df_map['TECAN_sample_target_position']
    return df_map
        
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
        
def write_pcr_report(df_map, outFH, mm_volume=0.0,
                     prm_volume=0.0, error_perc=10.0):
    """
    Writing a report on pcr reagnets
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
    outFH.write('\n# Volumes per RXN (ul)\n')
    write_report_line(outFH, 'Master Mix', mm_volume)
    write_report_line(outFH, 'Primers', prm_volume)
    ## raw total volumes
    outFH.write('\n# Total volumes (ul)\n')
    write_report_line(outFH, 'Master Mix', total_mm_volume)
    ## total volumes with error
    outFH.write('\n# Total volumes with (ul; {}% error)\n'.format(error_perc))
    write_report_line(outFH, 'Master Mix', total_mm_volume, error_perc=error_perc)
    # samples
    outFH.write('')


# main
if __name__ == '__main__':
    pass


