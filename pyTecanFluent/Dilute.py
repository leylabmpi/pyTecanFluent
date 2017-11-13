from __future__ import print_function
# import
## batteries
import os
import sys
import argparse
import functools
from itertools import cycle,product
## 3rd party
import numpy as np
import pandas as pd
## package
from pyTecanFluent import Utils
from pyTecanFluent import Fluent
from pyTecanFluent import Labware

# functions
def get_desc():
    desc = 'Create robot commands for diluting samples'
    return desc

def parse_args(test_args=None, subparsers=None):
    # desc
    desc = get_desc()
    epi = """DESCRIPTION:
    Create a worklist file for the TECAN Fluent robot for diluting samples.
    The input is an Excel or tab-delimited file the following columns:
    * "TECAN_labware_name" = Any name you want to give to your plate of samples
    * "TECAN_labware_type" = The labware type matching your samples
    * "TECAN_target_position" = The location of your samples in your labware (plate)
    * "TECAN_sample_conc" = The sample concentrations (numeric value; units=ng/ul) 
    
    Notes:
    * Sample locations in plates numbered are column-wise. 
    * All volumes are in ul.
    """
    if subparsers:
        parser = subparsers.add_parser('dilute', description=desc, epilog=epi,
                                       formatter_class=argparse.RawTextHelpFormatter)
    else:
        parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                         formatter_class=argparse.RawTextHelpFormatter)

    # args
    ## I/O
    groupIO = parser.add_argument_group('I/O')
    groupIO.add_argument('concfile', metavar='ConcFile', type=str,
                         help='An excel or tab-delim file of concentrations')
    groupIO.add_argument('--prefix', type=str, default='TECAN_dilute',
                         help='Output file name prefix (default: %(default)s)')

    ## concentration file
    conc = parser.add_argument_group('Concentation file')
    conc.add_argument('--format', type=str, default=None,
                        help='File format (excel or tab). If not provided, the format is determined from the file extension') 
    conc.add_argument('--header', action='store_false', default=True,
                        help='Header in the file? (default: %(default)s)')
    conc.add_argument('--rows', type=str, default='all',
                      help='Which rows (not including header) of the column file to use ("all"=all rows; "1-48"=rows 1-48) (default: %(default)s)')

    ## dilution
    dil = parser.add_argument_group('Dilution')
    dil.add_argument('--dilution', type=float, default=5.0,
                     help='Target dilution concentration (ng/ul) (default: %(default)s)')
    dil.add_argument('--minvolume', type=float, default=2.0,
                     help='Minimum sample volume to use (default: %(default)s)')
    dil.add_argument('--maxvolume', type=float, default=30.0,
                     help='Maximum sample volume to use (default: %(default)s)')
    dil.add_argument('--mintotal', type=float, default=10.0,
                     help='Minimum post-dilution total volume (default: %(default)s)')
    dil.add_argument('--dlabware_name', type=str, default='100ml_1',
                     help='Name of labware containing the dilutant (default: %(default)s)')
    dil.add_argument('--dlabware_type', type=str, default='100ml_1',
                     choices=['100ml_1', '1.5ml Eppendorf',
                              '2.0ml Eppendorf', '96 Well Eppendorf TwinTec PCR'], 
                     help='Labware type containing the dilutant (default: %(default)s)')

    ## destination plate
    dest = parser.add_argument_group('Destination labware')
    dest.add_argument('--destname', type=str, default='Diluted DNA plate',
                      help='Destination labware name (default: %(default)s)')
    dest.add_argument('--desttype', type=str, default='96 Well Eppendorf TwinTec PCR',
                      choices=['96 Well Eppendorf TwinTec PCR', '384 Well Biorad PCR'],                          
                      help='Destination labware type  on TECAN worktable (default: %(default)s)')

    ## tip type
    tips = parser.add_argument_group('Tip type')
    tips.add_argument('--tip1000_type', type=str, default='FCA, 1000ul SBS',
                      help='1000ul tip type (default: %(default)s)')
    tips.add_argument('--tip200_type', type=str, default='FCA, 200ul SBS',
                      help='200ul tip type (default: %(default)s)')
    tips.add_argument('--tip50_type', type=str, default='FCA, 50ul SBS',
                      help='50ul tip type (default: %(default)s)')
    tips.add_argument('--tip10_type', type=str, default='FCA, 10ul SBS',
                      help='10ul tip type (default: %(default)s)')
        
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
    df_conc = conc2df(args.concfile, 
                      file_format=args.format,
                      row_select=args.rows, 
                      header=args.header)

    # gwl object init 
    TipTypes = TipTypes={1000 : args.tip1000_type,
                         200 : args.tip200_type,
                         50 : args.tip50_type,
                         10 : args.tip10_type}    
    gwl = Fluent.gwl(TipTypes)
    
    # Determining dilution volumes
    df_conc = dilution_volumes(df_conc, 
                               dilute_conc=args.dilution,
                               min_vol=args.minvolume,
                               max_vol=args.maxvolume,
                               min_total=args.mintotal,
                               dest_type=args.desttype)
    
    # Adding destination data
    gwl.db.get_labware(args.desttype)
    n_wells = gwl.db.get_labware_wells(args.desttype)
    df_conc = add_dest(df_conc,
                       dest_name=args.destname,
                       dest_type=args.desttype,
                       dest_wells=n_wells)
            

    # Reordering dest if plate type is 384-well
    if n_wells == '384':
        df_conc = reorder_384well(df_conc, 'TECAN_dest_target_position')

    ## Dilutant
    pip_dilutant(df_conc, gwl=gwl,
                 src_labware_name=args.dlabware_name,
                 src_labware_type=args.dlabware_type)
    ## Sample
    pip_samples(df_conc, gwl=gwl)
    
    ## writing out worklist (gwl) file
    gwl_file = args.prefix + '.gwl'
    gwl.write(gwl_file)
    
    # making labware table
    lw = Labware.labware()
    lw.add_gwl(gwl)
    lw_df = lw.table()
    lw_file = args.prefix + '_labware.txt'
    lw_df.to_csv(lw_file, sep='\t', index=False)
    
    # Writing conc. out table
    conc_file = args.prefix + '_conc.txt'
    df_conc.round(1).to_csv(conc_file, sep='\t', index=False)

    # status
    Utils.file_written(gwl_file)
    Utils.file_written(conc_file)
    Utils.file_written(lw_file)

    # end
    return gwl_file, conc_file, lw_file


def check_args(args):
    """Checking user input
    """
    # input table column IDs
    args.rows = Utils.make_range(args.rows, set_zero_index=True)
    # dilution
    assert args.dilution >= 0.0, '--dilution must be >= 0'
    assert args.minvolume >= 0.0, '--minvolume must be >= 0'
    assert args.maxvolume > 0.0, '--maxvolume must be > 0'
    # tip type
    if args.tip1000_type.lower() == 'none':
        args.tip1000_type = None
    if args.tip200_type.lower() == 'none':
        args.tip200_type = None
    if args.tip50_type.lower() == 'none':
        args.tip50_type = None
    if args.tip10_type.lower() == 'none':
        args.tip10_type = None
        
def conc2df(concfile, row_select=None, file_format=None, header=True):
    """Loading a concentration file as a pandas dataframe
    """
    if header==True:
        header=0
    else:
        header=None
    # format
    if file_format is None:
        if concfile.endswith('.csv'):
            file_format = 'csv'
        elif concfile.endswith('.txt'):
            file_format = 'tab'
        elif concfile.endswith('.xls') or concfile.endswith('.xlsx'):
            file_format = 'excel'
    else:
        file_format = file_format.lower()
        
    # load via pandas IO
    if file_format == 'csv':
        df = pd.read_csv(concfile, sep=',', header=header)        
    elif file_format == 'tab':
        df = pd.read_csv(concfile, sep='\t', header=header)
    elif file_format == 'excel':
        xls = pd.ExcelFile(concfile)
        df = pd.read_excel(xls, header=header)
    else:
        raise ValueError('Concentration file not in usable format')

    # selecting rows
    if row_select is not None:
        df = df.iloc[row_select]
    
    # checking file format
    check_df_conc(df)

    # return
    return df

def missing_cols(df, req_cols):
    msg = 'Required column "{}" not found'
    for req_col in req_cols:
        if req_col not in df.columns.values:
            raise ValueError(msg.format(req_col))    

def check_df_conc(df_conc):
    """Assertions of df_conc object formatting
    """
    # checking for columns
    req_cols = ['TECAN_labware_name', 'TECAN_labware_type',
                'TECAN_target_position', 'TECAN_sample_conc']
    missing_cols(df_conc, req_cols)

    # checking labware types
    db = Fluent.db()
    msg = 'ERROR (concfile, line={}): labware type not recognized: {}'
    for i,lt in enumerate(df_conc['TECAN_labware_type']):
        try:
            db.get_labware(lt)
        except KeyError:
            raise KeyError(msg.format(i, lt))
                         
    # checking sample locations (>=1)
    msg = 'ERROR (concfile, line={}): location is < 1'
    for i,loc in enumerate(df_conc['TECAN_target_position']):
        if loc < 1:
            print(msg.format(i), file=sys.stderr)
    
    # checking sample conc
    msg = 'WARNING (concfile, line={}): concentration is <= 0'
    for i,sc in enumerate(df_conc['TECAN_sample_conc']):
        if sc <= 0.0:
            print(msg.format(i), file=sys.stderr)
            
def calc_sample_volume(row, dilute_conc, min_vol, max_vol):
    """sample_volume = dilute_conc * total_volume / conc 
    (v1 = c2*v2/c1)
    If sample_volume > max possibl volume to use, then just use max
    """
    # return 0 if conc <= 0 (this will be skipped)
    if row['TECAN_sample_conc'] <= 0:
        return 0
    # calc volume to use
    x = dilute_conc * row['TECAN_total_volume'] / row['TECAN_sample_conc']
    # ceiling
    if x > max_vol:
        x = max_vol
    # floor
    if x < min_vol:
        x = min_vol
    return x

def calc_dilutant_volume(row):
    """ dilutatant volume = total_volume - sample_volume
    """
    if row['TECAN_sample_volume'] <= 0:
        return 0
    x = row['TECAN_total_volume'] - row['TECAN_sample_volume']
    if x < 0:
        x = 0
    return x

def calc_total_volume(row, min_vol, max_vol, dilute_conc):
    """Calculating post-dlution volume
    """
    x = row['TECAN_sample_conc'] * min_vol / dilute_conc
    if x > max_vol:
        x = max_vol
    return x    

def calc_final_conc(row):
    """Calculating final conc (post-dilution sample conc.
    """
    if row['TECAN_sample_volume'] <= 0:
        return 0
    x = row['TECAN_sample_conc'] * row['TECAN_sample_volume']
    x = x / row['TECAN_total_volume']
    return x
    
def dilution_volumes(df_conc, dilute_conc, min_vol, max_vol, 
                     min_total, dest_type):
    """Setting the amoutn of sample to aliquot for dilution
    df_conc: pd.dataframe
    dilute_conc: concentration to dilute to 
    min_vol: min volume of sample to use
    max_vol: max total volume to use
    min_total: minimum total post-dilution volume
    dest_type: labware type for destination labware
    """
    # max well volume
    db = Fluent.db()
    db.get_labware(dest_type)
    max_well_vol = db.get_labware_max_volume(dest_type)

    # converting all negative concentrations to zero
    f = lambda row: 0 if row['TECAN_sample_conc'] < 0 else row['TECAN_sample_conc']
    df_conc['TECAN_sample_conc'] = df_conc.apply(f, axis=1)
    
    # range of dilutions
    samp_vol_range = max_vol - min_vol
    target_total_vol = round(samp_vol_range / 2 + min_vol)
    
    # final volume
    f = functools.partial(calc_total_volume, dilute_conc=dilute_conc,
                          min_vol=min_vol, max_vol=max_vol)
    df_conc['TECAN_total_volume'] = df_conc.apply(f, axis=1)
    if max(df_conc['TECAN_total_volume']) > max_well_vol:
        msg = 'ERROR: post-dilution volume exceeds max possible well volume.'
        msg += ' Lower --minvolume or chane destination labware type.'
        raise ValueError(msg)
    
    # raising total post-dilute volume if too low of dilute volume (if small dilution factor)
    df_conc.loc[df_conc.TECAN_total_volume < min_total, 'TECAN_total_volume'] = min_total
    # setting volumes
    f = functools.partial(calc_sample_volume, dilute_conc=dilute_conc,
                          min_vol=min_vol, max_vol=max_vol)
    df_conc['TECAN_sample_volume'] = df_conc.apply(f, axis=1)
    # dilutatant volume = total_volume - sample_volume
    df_conc['TECAN_dilutant_volume'] = df_conc.apply(calc_dilutant_volume, axis=1)
    # updating total volume
    f = lambda row: row['TECAN_sample_volume'] + row['TECAN_dilutant_volume']
    df_conc['TECAN_total_volume'] = df_conc.apply(f, axis=1)
    # calculating final conc
    df_conc['TECAN_final_conc'] = df_conc.apply(calc_final_conc, axis=1)
    ## target conc hit?
    msg_low = 'WARNING: (concfile, line{}): final concentration is low: {}'
    msg_high = 'WARNING: (concfile, line{}): final concentration is high: {}'
    for i,fc in enumerate(df_conc['TECAN_final_conc']):
        fc = round(fc, 1)
        if fc < round(dilute_conc, 1):            
            print(msg_low.format(i, fc), file=sys.stderr)
        if fc > round(dilute_conc, 1):
            print(msg_high.format(i, fc), file=sys.stderr)

    # return
    return df_conc
        
def add_dest(df_conc, dest_name, dest_type,  dest_wells=96):
    """Setting destination locations for samples & primers.
    Adding to df_conc:
      [dest_labware, dest_location]
    """
    dest_wells = int(dest_wells)

    # creating dest_name column; possibly multiple names
    n_dest_plates = int(round(df_conc.shape[0] / dest_wells + 0.5, 0))
    
    ## destination plate names
    if n_dest_plates > 1:
        dest_names = []
        for i in range(df_conc.shape[0]):
            x = int(round(i / dest_wells + 0.50001, 0))
            dest_names.append(dest_name + '[{:0>3}]'.format(x))   
        df_conc['TECAN_dest_labware_name'] = dest_names
    else:
        df_conc['TECAN_dest_labware_name'] = dest_name            

    ## labware type
    df_conc['TECAN_dest_labware_type'] = dest_type
    ## positions
    positions = cycle([x+1 for x in range(dest_wells)])
    positions = [next(positions) for x in range(df_conc.shape[0])]
    df_conc['TECAN_dest_target_position'] = positions

    # return
    return df_conc

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

def pip_dilutant(df_conc, gwl, src_labware_name, 
                 src_labware_type=None, lw_tracker=None):
    """Writing worklist commands for aliquoting dilutant.
    Using 1-asp-multi-disp with 200 ul tips.
    Method:
    * calc max multi-dispense for 50 or 200 ul tips 
    """
    # determing how many multi-disp per tip
    max_vol = max(df_conc.TECAN_dilutant_volume)
    tip_frac = 0.80 
    if max_vol > 1000 * tip_frac:
        raise ValueError('Max dilutant volume >{}ul'.format(1000 & tip_frac))
    if gwl.TipType_exists(10) and max_vol * 2 < 10 * tip_frac:
        n_disp = int(np.floor(10 * tip_frac / max_vol))   # using 10 ul tips
    elif gwl.TipType_exists(50) and max_vol * 2 < 50 * tip_frac:
        n_disp = int(np.floor(50 * tip_frac / max_vol))   # using 50 ul tips
    elif gwl.TipType_exists(200) and max_vol * 2 < 200 * tip_frac:
        n_disp = int(np.floor(200 * tip_frac / max_vol))  # using 200 ul tips
    elif gwl.TipType_exists(1000):
        n_disp = int(np.floor(1000 * tip_frac / max_vol))  # using 1000 ul tips
    else:
        raise ValueError('No TipType available for volume: {}ul'.format(max_vol))
        
    # making multi-disp object
    gwl.add(Fluent.comment('Dilutant'))
    MD = Fluent.multi_disp()
    MD.SrcRackLabel = src_labware_name
    MD.SrcRackType = src_labware_type
    MD.SrcPosition = 1                                   
    MD.DestRackLabel = df_conc.TECAN_dest_labware_name
    MD.DestRackType = df_conc.TECAN_dest_labware_type
    MD.DestPositions = df_conc.TECAN_dest_target_position
    MD.Volume = df_conc.TECAN_dilutant_volume             
    MD.NoOfMultiDisp = n_disp
    MD.add(gwl, tip_frac)
    
def pip_samples(df_conc, gwl):
    """Commands for aliquoting samples into dilutant
    """
    gwl.add(Fluent.comment('Samples'))
    # for each Sample-PCR_rxn_rep, write out asp/dispense commands
    for i in range(df_conc.shape[0]):
        # skipping no-volume
        if df_conc.loc[i,'TECAN_sample_volume'] <= 0:
            continue
        
        # aspiration
        asp = Fluent.aspirate()
        asp.RackLabel = df_conc.loc[i,'TECAN_labware_name']
        asp.RackType = df_conc.loc[i,'TECAN_labware_type']
        asp.Position = df_conc.loc[i,'TECAN_target_position']
        asp.Volume = round(df_conc.loc[i,'TECAN_sample_volume'], 2)
        asp.LiquidClass = 'Water Free Single No-cLLD'
        gwl.add(asp)
        
        # dispensing
        disp = Fluent.dispense()
        disp.RackLabel = df_conc.loc[i,'TECAN_dest_labware_name']
        disp.RackType = df_conc.loc[i,'TECAN_dest_labware_type']        
        disp.Position = df_conc.loc[i,'TECAN_dest_target_position']
        disp.Volume = round(df_conc.loc[i,'TECAN_sample_volume'], 2)
        disp.LiquidClass = 'Water Free Single No-cLLD'
        gwl.add(disp)

        # waste
        gwl.add(Fluent.waste())

# main
if __name__ == '__main__':
    pass


