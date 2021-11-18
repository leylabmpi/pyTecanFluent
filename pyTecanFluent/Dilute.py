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
    dil.add_argument('--min-volume', type=float, default=2.0,
                     help='Minimum sample volume to use (default: %(default)s)')
    dil.add_argument('--max-volume', type=float, default=30.0,
                     help='Maximum sample volume to use (default: %(default)s)')
    dil.add_argument('--min-total', type=float, default=10.0,
                     help='Minimum post-dilution total volume (default: %(default)s)')
    dil.add_argument('--only-dil', action='store_true', default=False,
                      help='If sample conc. is <=0, only add dilutant (default: %(default)s)')    
    dil.add_argument('--dil-labware-name', type=str, default='Dilutant',
                     help='Name of labware containing the dilutant (default: %(default)s)')
    dil.add_argument('--dil-labware-type', type=str, default='25ml_1 waste',
                     help='Labware type containing the dilutant (default: %(default)s)')
    dil.add_argument('--dil-liq', type=str, default='Water Free Single Wall Disp',
                      help='Dilutant liquid class (default: %(default)s)')
    dil.add_argument('--samp-liq', type=str, default='Water Free Single Wall Disp Aspirate Anyway',
                      help='Sample liquid class (default: %(default)s)')
    dil.add_argument('--reuse-tips', action='store_true', default=False,
                      help='Re-use tips for each dispense of dilutant? (default: %(default)s)')
        
    ## destination plate
    dest = parser.add_argument_group('Destination labware')
    dest.add_argument('--dest-name', type=str, default='Diluted sample plate',
                      help='Destination labware name (default: %(default)s)')
    dest.add_argument('--dest-type', type=str, default='96 Well Eppendorf TwinTec PCR',
                      help='Destination labware type on TECAN worktable (default: %(default)s)')
        
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
    TipTypes = ['FCA, 1000ul SBS', 'FCA, 200ul SBS',
                'FCA, 50ul SBS', 'FCA, 10ul SBS']    
    gwl = Fluent.gwl(TipTypes)
    
    # Determining dilution volumes
    df_conc = dilution_volumes(df_conc, 
                               dilute_conc=args.dilution,
                               min_vol=args.min_volume,
                               max_vol=args.max_volume,
                               min_total=args.min_total,
                               dest_type=args.dest_type,
                               only_dilutant=args.only_dil)
   
    # Adding destination data
    gwl.db.get_labware(args.dest_type)
    n_wells = gwl.db.get_labware_wells(args.dest_type)
    df_conc = add_dest(df_conc,
                       dest_name=args.dest_name,
                       dest_type=args.dest_type,
                       dest_wells=n_wells)
    df_conc = check_rack_labels(df_conc)
    
    # Reordering dest if plate type is 384-well
    df_conc = Utils.reorder_384well(df_conc, gwl,
                                   labware_name_col='TECAN_dest_labware_name',
                                   labware_type_col='TECAN_dest_labware_type',
                                   position_col='TECAN_dest_target_position')
        
    ## Dilutant
    pip_dilutant(df_conc, gwl=gwl,
                 src_labware_name=args.dil_labware_name,
                 src_labware_type=args.dil_labware_type,
                 liq_cls=args.dil_liq,
                 reuse_tips=args.reuse_tips)
    ## Sample
    pip_samples(df_conc, gwl=gwl,
                liq_cls=args.samp_liq)
    
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
    df_conc.round(2).to_csv(conc_file, sep='\t', index=False)

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
    assert args.min_volume >= 0.0, '--min-volume must be >= 0'
    assert args.max_volume > 0.0, '--max-volume must be > 0'
    # removing "tube" from end of labware type (if present)
    args.dil_labware_type = Utils.rm_tube(args.dil_labware_type)
    # special characters for namings
    args.dil_labware_name = Utils.rm_special_chars(args.dil_labware_name)
    args.dest_name = Utils.rm_special_chars(args.dest_name)
        
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
        elif concfile.endswith('.tsv'):
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

    # making sure labware names are "TECAN worklist friendly"
    df = Utils.rm_special_chars(df, 'TECAN_labware_name')

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
    
def check_rack_labels(df_conc):
    """Removing '.' for rack labels (causes execution failures)
    """
    cols = ['TECAN_labware_name', 'TECAN_dest_labware_name']
    for x in cols:
        df_conc[x] = [y.replace('.', '_') for y in df_conc[x].tolist()]
    return df_conc
            
def calc_final_volume(row, dilute_conc, min_vol, max_vol, min_total, max_total): 
    """Calculating post-dlution volume.
    row : row in pd.dataframe  (c1 = row['TECAN_sample_conc'])
    dilute_conc: target concentration to dilute to
    min_vol: min volume of sample to use  (min v1)
    max_vol: max volume of sample to use  (max v1)
    min_total: minimum total post-dilution volume (min v2)
    max_total: maximum total post-dilution volume (max v2)

    # algorithm (v2 = ?)
    # v2 = (v1 * c1) / c2
    ## min sample required to reach >= min_total while <= max_total
    """
    c1 = row['TECAN_sample_conc']
    # if can't dilute  take just sample (min_total volume)
    if c1 <= dilute_conc:
        msg = 'WARNING: (concfile, line{}): cannot dilute due to low conc.'
        print(msg.format(row.name + 1), file=sys.stderr)
        row['TECAN_total_volume'] = min_total
        row['TECAN_sample_volume'] = min_total
        row['TECAN_final_conc'] = c1
        return row
    # stepping from min to max possible sample volumes 
    v2_all = []
    c2_all = []
    for v1 in np.linspace(min_vol, max_vol, int((max_vol - min_vol) * 10)):
        v1 = round(v1, 2)
        v2 = round((v1 * c1) / dilute_conc, 2)
        if v2 == 0:
            c2 = 0
        else:
            c2 = round((v1 * c1) / v2, 2)
        if v2 >= min_total and v2 <= max_total:
            row['TECAN_total_volume'] = v2
            row['TECAN_sample_volume'] = v1
            row['TECAN_final_conc'] = c2
            return row
        else:
            v2_all.append(v2)
            c2_all.append(c2)
    # if no volumes work, then use volume with closest to conc
    ## Note: must be >= min_total
    min_delta_c = [abs(dilute_conc - c2) for c2 in c2_all]
    c2 = c2_all[min_delta_c.index(min(min_delta_c))]
    v2 = v2_all[min_delta_c.index(min(min_delta_c))]
    v2 = v2 if v2 >= min_total else min_total
    v2 = round(v2, 2)
    c2 = round((v1 * c1) / v2, 2)
    if c2 != dilute_conc:
        msg = 'WARNING: (concfile, line{}): final concentration is {}'
        print(msg.format(row.name + 1, c2), file=sys.stderr)
    row['TECAN_total_volume'] = v2
    row['TECAN_sample_volume'] = v1
    row['TECAN_final_conc'] = c2
    return row

def calc_dilutant_volume(row):
    """ dilutatant volume = total_volume - sample_volume
    """
    x =  row['TECAN_total_volume'] - row['TECAN_sample_volume']
    if x < 0:
        msg = 'ERROR: (concfile, line{}): dilutant volume = {}'
        raise ValueError(msg.format(row.name + 1, x))
    return x
    
def dilution_volumes(df_conc, dilute_conc, min_vol, max_vol, 
                     min_total, dest_type, only_dilutant=False):
    """Setting the amount of sample to aliquot for dilution
    df_conc: pd.dataframe
    dilute_conc: concentration to dilute to 
    min_vol: min volume of sample to use
    max_vol: max volume of sample to use
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
    
    # calc final volume
    f = functools.partial(calc_final_volume,
                          dilute_conc=dilute_conc,
                          min_vol=min_vol,
                          max_vol=max_vol,
                          min_total=min_total,
                          max_total=max_well_vol)
    df_conc = df_conc.apply(f, axis=1)
    if max(df_conc['TECAN_total_volume']) > max_well_vol:
        msg = 'ERROR: post-dilution volume exceeds max possible well volume.'
        msg += ' Lower --minvolume or chane destination labware type.'
        raise ValueError(msg)
        
    # calc dilutant volume (final_volume - sample_volume)
    f = functools.partial(calc_dilutant_volume)
    df_conc['TECAN_dilutant_volume'] = df_conc.apply(f, axis=1)
        
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

    # reorder columns
    cols = ['TECAN_labware_name',
            'TECAN_labware_type',
            'TECAN_target_position',
            'TECAN_dest_labware_name',
            'TECAN_dest_labware_type',
            'TECAN_dest_target_position',
            'TECAN_sample_conc',
            'TECAN_sample_volume',
            'TECAN_dilutant_volume',
            'TECAN_total_volume',
            'TECAN_final_conc']
    df_conc = df_conc[cols]
    
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

def pip_dilutant(df_conc, gwl, src_labware_name, src_labware_type=None,
                 liq_cls='Water Free Single', reuse_tips=False):
    """Commands for aliquoting dilutant.
    """
    gwl.add(Fluent.Comment('Dilutant'))
    # filtering by volume range
    DTH_vols = list(gwl.get_DTH_volumes().values())
    DTH_vols = sorted(DTH_vols)
    DTH_vols = zip([-1] + DTH_vols[:-1], DTH_vols)

    # filtering to just values for each DTH volume range
    for DTH_vol in DTH_vols:
        x = (df_conc['TECAN_dilutant_volume'] > DTH_vol[0]) & \
            (df_conc['TECAN_dilutant_volume'] <= DTH_vol[1])
        df_conc_tmp = df_conc.loc[x]
        if df_conc_tmp.shape[0] <= 0:
            continue
        df_conc_tmp.reset_index(inplace=True)
        # for each sample, transfer aliquot via asp-disp
        for i in range(df_conc_tmp.shape[0]):
            # skipping no-volume
            volume = round(df_conc_tmp.loc[i,'TECAN_dilutant_volume'], 2)
            if volume <= 0.0:
                continue        
            # aspiration
            asp = Fluent.Aspirate()
            asp.RackLabel = src_labware_name
            asp.RackType = src_labware_type
            asp.Position = 1
            asp.Volume = volume
            asp.LiquidClass = liq_cls
            gwl.add(asp)        
            # dispensing
            disp = Fluent.Dispense()
            disp.RackLabel = df_conc_tmp.loc[i,'TECAN_dest_labware_name']
            disp.RackType = df_conc_tmp.loc[i,'TECAN_dest_labware_type']        
            disp.Position = df_conc_tmp.loc[i,'TECAN_dest_target_position']
            disp.Volume = volume
            disp.LiquidClass = liq_cls
            gwl.add(disp)
            # end
            if reuse_tips is True:
                gwl.add(Fluent.Flush())
            else:
                gwl.add(Fluent.Waste())
        # waste
        if reuse_tips is True:
            gwl.add(Fluent.Waste())        
    # adding break
    gwl.add(Fluent.Break())

def pip_samples(df_conc, gwl, liq_cls='Water Free Single'):
    """Commands for aliquoting samples into dilutant
    """
    gwl.add(Fluent.Comment('Samples'))
    
    # for each sample, transfer aliquot via asp-disp 
    for i in range(df_conc.shape[0]):
        # skipping no-volume
        volume = round(df_conc.loc[i,'TECAN_sample_volume'], 2)
        if volume <= 0:
            continue        
        # aspiration
        asp = Fluent.Aspirate()
        asp.RackLabel = df_conc.loc[i,'TECAN_labware_name']
        asp.RackType = df_conc.loc[i,'TECAN_labware_type']
        asp.Position = df_conc.loc[i,'TECAN_target_position']
        asp.Volume = volume
        asp.LiquidClass = liq_cls
        gwl.add(asp)        
        # dispensing
        disp = Fluent.Dispense()
        disp.RackLabel = df_conc.loc[i,'TECAN_dest_labware_name']
        disp.RackType = df_conc.loc[i,'TECAN_dest_labware_type']        
        disp.Position = df_conc.loc[i,'TECAN_dest_target_position']
        disp.Volume = volume
        disp.LiquidClass = liq_cls
        gwl.add(disp)

        # waste
        gwl.add(Fluent.Waste())

# main
if __name__ == '__main__':
    pass


