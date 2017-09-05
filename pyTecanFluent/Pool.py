from __future__ import print_function
# import
## batteries
import os
import sys
import argparse
import functools
from itertools import product
## 3rd party
import numpy as np
import pandas as pd
## package
from pyTecanFluent import Utils
from pyTecanFluent import Fluent
from pyTecanFluent import Labware

# functions
def get_desc():
    desc = 'Create robot commands for pooling samples'
    return desc

def parse_args(test_args=None, subparsers=None):
    # desc
    desc = get_desc()
    epi = """DESCRIPTION:
    Create a worklist and labware file for the TECAN Fluent robot for pooling samples
    (eg., reaction replicates). 
    The input is an Excel or tab-delimited file containing 2 columns:
    1) A column of sample names (samples with the same name will be pooled)
    2) A column designating whether to include or skip the samplles 
       (include = 'Success/Pass/Include'; skip = 'Fail/Skip')
    
    Mapping file:
    If a mapping file is provided (same names as in the pooling file),
    then the mapping file will be trimmed to just those pooled, and
    the final pooled locations will be added to the mapping table. 

    Notes:
    * Sample locations in plates numbered are column-wise. 
    * All volumes are in ul.
    """
    if subparsers:
        parser = subparsers.add_parser('pool', description=desc, epilog=epi,
                                       formatter_class=argparse.RawTextHelpFormatter)
    else:
        parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                         formatter_class=argparse.RawTextHelpFormatter)


    #----------------------------------------------------------#
    # TODO: allow multiple sample files; need to provide plateID 
        
    # args
    ## I/O
    groupIO = parser.add_argument_group('I/O')
    groupIO.add_argument('samplefile', metavar='SampleFile', type=str,
                         help='An excel or tab-delim file of samples to pool')
    groupIO.add_argument('--prefix', type=str, default='TECAN_pool',
                         help='Output file name prefix (default: %(default)s)')
    groupIO.add_argument('--mapfile', type=str, 
                         help='A QIIME-formatted mapping file')
    
    ## file format
    filef = parser.add_argument_group('sample file')
    ### sample file
    filef.add_argument('--sample_format', type=str, default=None,
                        help='File format (excel or tab). If not provided, the format is determined from the file extension') 
    filef.add_argument('--sample_header', action='store_false', default=True,
                        help='Header in the sample file? (default: %(default)s)')
    filef.add_argument('--sample_rows', type=str, default='all',
                      help='Which rows (not including header) of the sample file to use ("all"=all rows; "1-48"=rows 1-48) (default: %(default)s)')
    filef.add_argument('--sample_col', type=str, default='Sample',
                        help='Column containing the samples  (default: %(default)s)')
    filef.add_argument('--include_col', type=str, default='Call',
                        help='Column designating sample include/skip? (default: %(default)s)')
    filef.add_argument('--pos_col', type=str, defaul='Well',
                        help='Column designating sample location in the plate')
    ### mapping file
    filef.add_argument('--map_format', type=str, default=None,
                        help='File format (excel or tab). If not provided, the format is determined from the file extension') 
    filef.add_argument('--map_header', action='store_false', default=True,
                        help='Header in the mapping file? (default: %(default)s)')

    ## pooling
    pooling = parser.add_argument_group('Pooling')
    pooling.add_argument('--volume', type=float, default=30.0,
                         help='Per-sample volume to pool (default: %(default)s)')
    pooling.add_argument('--liqcls', type=str, default='Water Contact Free Multi No-cLLD',
                         help='Liquid class for pooling (default: %(default)s)')

    ## destination plate
    dest = parser.add_argument_group('Destination labware')
    dest.add_argument('--destname', type=str, default='Diluted DNA plate',
                      help='Destination labware name (default: %(default)s)')
    dest.add_argument('--desttype', type=str, default='96 Well Eppendorf TwinTec PCR',
                      choices=['96 Well Eppendorf TwinTec PCR', '384 Well Biorad PCR'],                          
                      help='Destination labware type  on TECAN worktable (default: %(default)s)')
    dest.add_argument('--deststart', type=int, default=1,
                      help='Starting position (well) on the destination labware (default: %(default)s)')

    
    # parse & return
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
    ## sample file
    df_samp = sample2df(args.samplefile, 
                        sample_col=args.sample_col,
                        include_col=args.include_col,
                        file_format=args.sample_format,
                        row_select=args.sample_rows, 
                        header=args.sample_header)

    ## mapping file
    df_map = map2df(args.mapfile,
                    file_format=args.map_format,
                    header=args.map_header)

    # filtering sample file to just pass
    df_samp = df_samp.loc[df_samp[args.include_col].isin(['success', 'pass', 'include'])]
    print(df_samp)
    
    # # Reordering dest if plate type is 384-well
    # try:
    #     n_wells = Labware.LABWARE_DB[args.desttype]['wells']
    # except KeyError:
    #     msg = 'Labware type "{}" does not have "wells" attribute'
    #     raise KeyError(msg.format(args.desttype))
    # if n_wells == '384':
    #     df_conc = reorder_384well(df_conc, 'TECAN_dest_target_position')

    # # Writing out gwl file
    # lw_tracker = Labware.labware_tracker()
    # gwl_file = args.prefix + '.gwl'
    # with open(gwl_file, 'w') as gwlFH:
    #     ## Dilutant
    #     pip_dilutant(df_conc, outFH=gwlFH,
    #                  src_labware_name=args.dlabware_name,
    #                  src_labware_type=args.dlabware_type,
    #                  lw_tracker=lw_tracker)
    #     ## Sample
    #     pip_samples(df_conc, outFH=gwlFH,
    #                 lw_tracker=lw_tracker)

    # # making labware table
    # df_labware = lw_tracker.labware_table()
    # lw_file = args.prefix + '_labware.txt'
    # df_labware.to_csv(lw_file, sep='\t', index=False)

    # # Writing out table
    # conc_file = args.prefix + '_conc.txt'
    # df_conc.round(1).to_csv(conc_file, sep='\t', index=False)

    # # Create windows-line breaks formatted versions
    # gwl_file_win = Utils.to_win(gwl_file)
    # conc_file_win = Utils.to_win(conc_file)
    # lw_file_win = Utils.to_win(lw_file)

    # # status
    # Utils.file_written(gwl_file)
    # Utils.file_written(conc_file)
    # Utils.file_written(lw_file)
    # Utils.file_written(gwl_file_win)
    # Utils.file_written(conc_file_win)    
    # Utils.file_written(lw_file_win)

    
    # # end
    # return (gwl_file, gwl_file_win, 
    #         conc_file, conc_file_win, 
    #         lw_file, lw_file_win)


def check_args(args):
    """Checking user input
    """
    # input table column IDs
    args.rows = Utils.make_range(args.sample_rows, set_zero_index=True)
    # dilution
    assert args.volume >= 0.0, '--volume must be >= 0'
    # destination labware type
    try:
        Labware.LABWARE_DB[args.desttype]
    except KeyError:
        msg = 'Destination labware type "{}" not recognized'
        raise ValueError(msg.format(args.desttype))

                         
def sample2df(samplefile, sample_col, include_col, row_select=None, file_format=None, header=True):
    """Loading a sample file as a pandas dataframe
    """
    if header==True:
        header=0
    else:
        header=None
    # format
    if file_format is None:
        if samplefile.endswith('.csv'):
            file_format = 'csv'
        elif samplefile.endswith('.txt'):
            file_format = 'tab'
        elif samplefile.endswith('.xls') or samplefile.endswith('.xlsx'):
            file_format = 'excel'
    else:
        file_format = file_format.lower()
        
    # load via pandas IO
    if file_format == 'csv':
        df = pd.read_csv(samplefile, sep=',', header=header)        
    elif file_format == 'tab':
        df = pd.read_csv(samplefile, sep='\t', header=header)
    elif file_format == 'excel':
        xls = pd.ExcelFile(samplefile)
        df = pd.read_excel(xls, header=header)
    else:
        raise ValueError('Sample file is not in a usable format')

    # checking dataframe format
    ## columns
    req_cols = [sample_col, include_col]
    msg = 'Column "{}" not found in sample table'
    for req_col in req_cols:
        if req_col not in df.columns.values:
            raise ValueError(msg.format(req_col))    
    ## include col
    f = lambda x: x.lower()
    df.ix[:,include_col] = df.ix[:,include_col].apply(f)
    msg = '"{}" value not allowed in include column in sample file'
    df.ix[:,include_col].apply(check_include_column)
    
    # selecting relevant columns
    df = df.ix[:,[sample_col, include_col]]

    
    # return
    return df

def check_include_column(x):
    psbl_vals = ('success', 'include', 'pass', 'fail', 'skip')
    assert x in psbl_vals, msg.format(x)

def map2df(mapfile, file_format=None, header=True):
    """Loading a mapping file as a pandas dataframe
    """
    if mapfile is None:
        return None
    if header==True:
        header=0
    else:
        header=None
    # format
    if file_format is None:
        if mapfile.endswith('.csv'):
            file_format = 'csv'
        elif mapfile.endswith('.txt'):
            file_format = 'tab'
        elif mapfile.endswith('.xls') or mapfile.endswith('xlsx'):
            file_format = 'excel'
    else:
        file_format = file_format.lower()
        
    # load via pandas IO
    if file_format == 'csv':
        df = pd.read_csv(mapfile, sep=',', header=header)        
    elif file_format == 'tab':
        df = pd.read_csv(mapfile, sep='\t', header=header)
    elif file_format == 'excel':
        xls = pd.ExcelFile(mapfile)
        df = pd.read_excel(xls, header=header)
    else:
        raise ValueError('Mapping file is not in a usable format')
    
    # checking for required columns
    req_cols = ['#SampleID']
    msg = 'Column "{}" not found in mapping file'
    for req_col in req_cols:
        if req_col not in df.columns.values:
            raise ValueError(msg.format(req_col))    
    
    # ret
    return df

    
# main
if __name__ == '__main__':
    pass


