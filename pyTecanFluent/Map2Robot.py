from __future__ import print_function
# import
## batteries
import os
import sys
import argparse
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
    desc = 'Convert a mapping file to a NGS amplicon worklist file for the TECAN robot'
    return desc

def parse_args(test_args=None, subparsers=None):
    # desc
    #desc = 'Convert a mapping file to a NGS amplicon worklist file for the TECAN robot'
    desc = get_desc()
    epi = """DESCRIPTION:
    Convert a QIIME-formatted mapping file to a GWL file, which is used by the TECAN
    robot to conduct the NGS amplicon PCR prep (ie., combining MasterMix, primers, samples, etc).
    The extra columns in the mapping file designate the SOURCE of samples and primers;
    the DESTINATION (plate & well) is set by this script. 

    EXTRA COLUMNS in MAPPING FILE:

    # If no barcoding, single barcoding, or dual barcodes already combined:
    * "TECAN_sample_labware_name" = The sample labware name on the robot worktable
    * "TECAN_sample_labware_type" = The type of labware containing samples (eg., '96 Well Eppendorf TwinTec PCR')
    * "TECAN_sample_target_position" = The well or tube location (a number)
    * "TECAN_sample_rxn_volume" = The volume of sample to use per PCR (ul)
    * "TECAN_primer_labware_name" = The primer plate labware name on the robot worktable
    * "TECAN_primer_labware_type" = The primer plate labware type on the robot worktable
    * "TECAN_primer_target_position" = The well location (1-96 or 1-384)

    CONTROLS:
    * For the positive & negative controls, include them in the mapping file.
    * If the controls (or samples) are provided in a tube, use "micro15[XXX]" 
      for the "TECAN_sample_labware" column, but change "XXX" to the tube number
      that you want to use (eg., micro15[003] for tube position 3)

    OUTPUT FILES:
    * The output files ending in "_win" have Windows line breads (needed for the robot)

    MISC NOTES:
    * All volumes are in ul
    * Plate well locations are 1 to n-wells; numbering by column
    * PicoGreen should be added to the MasterMix *prior* to loading on robot
    """
    if subparsers:
        parser = subparsers.add_parser('map2robot', description=desc, epilog=epi,
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
    groupIO.add_argument('--prefix', type=str, default='TECAN_NGS_amplicon',
                         help='Output file name prefix (default: %(default)s)')

    ## Destination plate
    dest = parser.add_argument_group('Destination plate')
    dest.add_argument('--destname', type=str, default='Destination plate',
                      help='Distination labware name (default: %(default)s)')
    dest.add_argument('--desttype', type=str, default='96 Well Eppendorf TwinTec PCR',
                      choices=['96 Well Eppendorf TwinTec PCR',
                               '384 Well Biorad PCR'],
                      help='Destination labware type (default: %(default)s)')
    dest.add_argument('--deststart', type=int, default=1,
                      help='Start well number on destination plate (default: %(default)s)')
    dest.add_argument('--rxns', type=int, default=3,
                      help='Number of replicate PCRs per sample (default: %(default)s)')

    ## MasterMix
    mm = parser.add_argument_group('Master mix')
    mm.add_argument('--mmtube', type=int, default=1,
                        help='MasterMix tube number (default: %(default)s)')
    mm.add_argument('--mmvolume', type=float, default=13.1,
                        help='MasterMix volume per PCR (default: %(default)s)')

    ## Primers
    primers = parser.add_argument_group('Primers')
    primers.add_argument('--fpvolume', type=float, default=1.0,
                        help='Forward primer volume per PCR (default: %(default)s)')
    primers.add_argument('--rpvolume', type=float, default=1.0,
                        help='Reverse primer volume per PCR (default: %(default)s)')
    primers.add_argument('--fptube', type=int, default=0,
                        help='Forward non-bacode primer tube number (0 = barcoded primer on a plate), (default: %(default)s)')
    primers.add_argument('--rptube', type=int, default=0,
                        help='Reverse non-bacode primer tube number (0 = barcoded primer on a plate), (default: %(default)s)')
    
    # Liquid classes
    liq = parser.add_argument_group('Liquid classes')
    liq.add_argument('--mm-liq', type=str, default='MasterMix Free Multi',
                      help='Mastermix liquid class (default: %(default)s)')
    liq.add_argument('--primer-liq', type=str, default='Water Contact Wet Single',
                      help='Primer liquid class (default: %(default)s)')
    liq.add_argument('--sample-liq', type=str, default='Water Contact Wet Single',
                      help='Sample liquid class (default: %(default)s)')
    liq.add_argument('--water-liq', type=str, default='Water Contact Wet Single',
                      help='Water liquid class (default: %(default)s)')

    ## Misc
    misc = parser.add_argument_group('Misc')
    misc.add_argument('--pcrvolume', type=float, default=25.0,
                        help='Total volume per PCR (default: %(default)s)')
    misc.add_argument('--errorperc', type=float, default=10.0,
                        help='Percent of extra total reagent volume to include (default: %(default)s)')

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
    df_map = map2df(args.mapfile, row_select=args.rows)
    check_df_map(df_map, args)

    # Making destination dataframe
    df_map = add_dest(df_map, args.destname,
                      dest_type=args.desttype,
                      dest_start=args.deststart,
                      rxn_reps=args.rxns)

    positions = get_wells(args.desttype)
    # Reordering dest if plate type is 384-well
    if positions == '384':
        df_map = reorder_384well(df_map, 'TECAN_dest_target_position')

    # GWL file construction
    lw_tracker = Labware.labware_tracker()
    ## gwl open
    gwl_file = args.prefix + '.gwl'
    with open(gwl_file, 'w') as gwlFH:
        ## mastermix
        pip_mastermix(df_map, gwlFH,
                      mmtube=args.mmtube,
                      mmvolume=args.mmvolume, 
                      liq_cls=args.mm_liq, 
                      lw_tracker=lw_tracker)
        ## primers
        pip_primers(df_map, gwlFH,
                    fp_volume=args.fpvolume,
                    rp_volume=args.rpvolume,
                    fp_tube=args.fptube,
                    rp_tube=args.rptube, 
                    liq_cls=args.primer_liq, 
                    lw_tracker=lw_tracker)
        ## samples
        pip_samples(df_map, gwlFH, 
                    liq_cls=args.sample_liq, 
                    lw_tracker=lw_tracker)
        ## water
        pip_water(df_map, gwlFH,
                  pcr_volume=args.pcrvolume,
                  mm_volume=args.mmvolume,
                  fp_volume=args.fpvolume,
                  rp_volume=args.rpvolume, 
                  liq_cls=args.water_liq, 
                  lw_tracker=lw_tracker)

    # Report (total volumes; sample truncation; samples)
    report_file = args.prefix + '.report'
    with open(report_file, 'w') as repFH:
        write_report(df_map, outFH=repFH,
                     pcr_volume=args.pcrvolume,
                     mm_volume=args.mmvolume,
                     fp_tube=args.fptube,
                     fp_volume=args.fpvolume,
                     rp_tube=args.rptube,
                     rp_volume=args.rpvolume,
                     n_rxn_reps=args.rxns,
                     error_perc=args.errorperc)

    # making labware table
    df_labware = lw_tracker.labware_table()
    lw_file = args.prefix + '_labware.txt'
    df_labware.to_csv(lw_file, sep='\t', index=False)

    # Mapping file with destinations
    df_file = args.prefix + '_map.txt'
    df_map['TECAN_water_rxn_volume'] = df_map['TECAN_water_rxn_volume'].round(2)
    df_map['TECAN_dest_target_position'] = df_map['TECAN_dest_target_position'].astype(int)
    df_map['TECAN_pcr_rxn_rep'] = df_map['TECAN_pcr_rxn_rep'].astype(int)
    df_map.to_csv(df_file, sep='\t', index=False, na_rep='NA')

    # Create windows-line breaks formatted versions
    gwl_file_win = Utils.to_win(gwl_file)
    report_file_win = Utils.to_win(report_file)
    df_file_win = Utils.to_win(df_file)
    lw_file_win = Utils.to_win(lw_file)

    # status on files written
    Utils.file_written(gwl_file)
    Utils.file_written(report_file)
    Utils.file_written(df_file)
    Utils.file_written(lw_file)
    Utils.file_written(gwl_file_win)
    Utils.file_written(report_file_win)
    Utils.file_written(df_file_win)
    Utils.file_written(lw_file_win)
    
    # Return
    return (gwl_file, gwl_file_win,
            report_file, report_file_win, 
            df_file, df_file_win, 
            lw_file, lw_file_win)


def get_wells(dest_type):
    try:
        Labware.LABWARE_DB[dest_type]
    except KeyError:
        msg = 'Labware type "{}" not recognised'
        raise KeyError(msg.format(dest_type))
    try:
        positions = Labware.LABWARE_DB[dest_type]['wells']
    except KeyError:
        msg = 'No "wells" key for labware type "{}"'
        raise KeyError(msg.format(dest_type))
    return positions
    

def check_args(args):
    """Checking user input
    """
    # destination start
    positions = get_wells(args.desttype)
    if args.deststart < 1 or args.deststart > positions:
        msg = 'Destination start well # must be in range: 1-{}'
        raise ValueError(msg.format(destlimit))

    # rows in mapping file
    args.rows = Utils.make_range(args.rows, set_zero_index=True)
    # tube number
    if args.mmtube < 1 or args.mmtube > 24:
        msg = '{} tube # must be in range: 1-24'
        raise ValueError(msg.format('MasterMix'))
    if args.fptube < 0 or args.fptube > 24:
        msg = '{} tube # must be in range: 1-24 (or 0 if no tube)'
        raise ValueError(msg.format('Forward primer'))
    if args.rptube < 0 or args.rptube > 24:
        msg = '{} tube # must be in range: 1-24 (or 0 if no tube)'
        raise ValueError(msg.format('Reverse primer'))
    # volumes
    ## mastermix
    if args.mmvolume > args.pcrvolume:
        msg = 'MasterMix volume > total PCR volume'
        raise ValueError(msg)
    if args.mmvolume / 2.0 > args.pcrvolume:
        msg = 'WARNING: MasterMix volume > half of PCR volume'
        print(msg, file=sys.stderr)


def map2df(mapfile, row_select=None):
    """Loading a mapping file as a pandas dataframe
    mapfile: string; mapping file path
    row_select: select particular rows of table in map file
    """
    # load via pandas IO
    if mapfile.endswith('.txt') or mapfile.endswith('.csv'):
        df = pd.read_csv(mapfile, sep='\t')
    elif mapfile.endswith('.xls') or mapfile.endswith('.xlsx'):
        xls = pd.ExcelFile(mapfile)
        df = pd.read_excel(xls)
    else:
        raise ValueError('Mapping file not in usable format')

    # selecting particular rows
    if row_select is not None:
        df = df.iloc[row_select]

    # return
    return df

def missing_cols(df, req_cols):
    msg = 'Required column "{}" not found'
    for req_col in req_cols:
        if req_col not in df.columns.values:
            raise ValueError(msg.format(req_col))    

def check_df_map(df_map, args):
    """Assertions of df_map object formatting
    * Assumes `sample` field = 1st column
    df_map: map dataframe
    """
    # checking columns
    ## universal
    req_cols = ['TECAN_sample_labware_name', 'TECAN_sample_labware_type',
                'TECAN_sample_target_position', 'TECAN_sample_rxn_volume', 
                'TECAN_primer_labware_name', 'TECAN_primer_labware_type']
    missing_cols(df_map, req_cols)

    # checking for unique samples
    if any(df_map.duplicated(df_map.columns.values[0])):
        msg = 'WARNING: duplicated sample values in the mapping file'
        print(msg, file=sys.stderr)

    # checking for unique barcode locations (NaN's filtered out)
    dups = df_map[req_cols].dropna().duplicated(keep=False)
    if any(dups):
        msg = 'WARNING: duplicated barcodes in the mapping file'
        print(msg, file=sys.stderr)
        
    # checking sample volumes
    msg = 'WARNING: sample volume > mastermix volume'
    for sv in df_map['TECAN_sample_rxn_volume']:
        if sv > args.mmvolume:
            print(msg, file=sys.stderr)
        if sv < 0:
            raise ValueError('Sample volume < 0')


def add_dest(df_map, dest_labware, dest_type='96 Well Eppendorf TwinTec PCR',
             dest_start=1, rxn_reps=3):
    """Setting destination locations for samples & primers.
    Making a new dataframe with:
      [sample, sample_rep, dest_labware, dest_location]
    * For each sample (i):
      * For each replicate (ii):
        * plate = destination plate type
        * well = i * (ii+1) + (ii+1) + start_offset
    Joining to df_map
    """
    dest_start= int(dest_start)
    # labware type found in DB?
    positions = get_wells(dest_type)

    # init destination df
    sample_col = df_map.columns[0]
    cols = [sample_col, 'TECAN_pcr_rxn_rep',
            'TECAN_dest_labware_name',
            'TECAN_dest_labware_type', 
            'TECAN_dest_target_position']    
    ncol = len(cols)
    nrow = df_map.shape[0] * rxn_reps        # number of rxns
    #if nrow > positions + 1 - dest_start:
    #    nrow = positions + 1 - dest_start
    #elif nrow > positions + 1 - dest_start:
    #    nrow = positions + 1 - dest_start
    df_dest = pd.DataFrame(np.nan, index=range(nrow), columns=cols)

    # number of destination plates required
    n_dest_plates = round(df_dest.shape[0] / positions + 0.5, 0)
    if n_dest_plates > 1:
        msg = ('WARNING: Not enough wells for the number of samples.' 
        ' Using multiple destination plates')
        print(msg, file=sys.stderr)
        
    
    # filling destination df
    orig_dest_labware = dest_labware
    for i,(sample,rep) in enumerate(product(df_map.ix[:,0], range(rxn_reps))):
        # dest location
        dest_position = i + dest_start
        dest_position = positions if dest_position % positions == 0 else dest_position % positions 

        # destination plate name
        if n_dest_plates > 1:
            x = round((i + dest_start) / positions + 0.5, 0)        
            dest_labware = '{} {}'.format(orig_dest_labware, int(x))

        # adding values DF
        df_dest.iloc[i] = [sample, rep+1, dest_labware, dest_type, dest_position]

    # df join (map + destination)
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


def pip_mastermix(df_map, outFH, mmvolume=13.1, mmtube=1, 
                  liq_cls='Water Free Single', lw_tracker=None):
    """Writing worklist commands for aliquoting mastermix.
    Using 1-asp-multi-disp with 200 ul tips.
    Method:
    * calc max multi-dispense for 200 ul tips & mmvolume
    * for 1:n_dispense
      * determine how many disp left for channel (every 8th)
      * if n_disp_left < n_disp: n_disp = n_disp_left
      * calc total volume: n_disp * mmvolume
    """
    mmtube = int(mmtube)
    outFH.write('C;MasterMix\n')

    MD = Fluent.multi_disp(lw_tracker)
    MD.SrcRackLabel = 'Mastermix tube[{0:0>3}]'.format(mmtube)
    MD.SrcRackType = '1.5ml Eppendorf'
    MD.SrcPosition = mmtube
    MD.DestRackLabel = df_map.ix[:,'TECAN_dest_labware_name']    
    MD.DestRackType = df_map.ix[:,'TECAN_dest_labware_type']    
    MD.DestPositions = df_map.ix[:,'TECAN_dest_target_position']
    MD.Volume = mmvolume
    MD.LiquidClass = liq_cls
    MD.NoOfMultiDisp = int(np.floor(160 / mmvolume))  # using 200 ul tips

    outFH.write(MD.cmd() + '\n')


def pip_nonbarcode_primer(df_map, outFH, volume, tube, liq_cls='Water Free Single',
                          lw_tracker=None):
    """Pipetting primers from tube.
    Assuming primer is aliquoted to all samples
    df_map : mapping file dataframe
    volume : volume to aliquot to each reaction
    tube : tube number (rackID assumed)
    liq_cls : liquid class
    """
    tube = int(tube)
    outFH.write('C;Non-barcoded primers\n')
    # for each Sample-PCR_rxn_rep, write out asp/dispense commands
    for i in range(df_map.shape[0]):
        # aspiration
        asp = Fluent.aspirate()
        asp.RackLabel = 'Primer tube[{0:0>3}]'.format(tube)
        asp.RackType = '1.5ml Eppendorf'
        asp.Position = tube
        asp.Volume = volume
        asp.LiquidClass = liq_cls
        outFH.write(asp.cmd() + '\n')

        # dispensing
        disp = Fluent.dispense()
        disp.RackLabel = df_map.ix[i,'TECAN_dest_labware_name']
        disp.RackType = df_map.ix[i,'TECAN_dest_labware_type']
        disp.Position = df_map.ix[i,'TECAN_dest_target_position']
        disp.Volume = volume
        disp.LiquidClass = liq_cls
        outFH.write(disp.cmd() + '\n')

        # tip to waste
        outFH.write('W;\n')

        # tracking labware & tip usage
        lw_tracker.add(asp)
        lw_tracker.add(disp, add_tip=False)


def pip_primer(i, outFH, df_map, primer_labware_name, 
               primer_labware_type, primer_target_position,
               primer_plate_volume, liq_cls, lw_tracker=None):
    """Pipetting a single primer
    """
    # aspiration
    asp = Fluent.aspirate()
    asp.RackLabel = df_map.ix[i, primer_labware_name]
    asp.RackType = df_map.ix[i, primer_labware_type]
    asp.Position = df_map.ix[i, primer_target_position]
    asp.Volume = primer_plate_volume
    asp.LiquidClass = liq_cls
    outFH.write(asp.cmd() + '\n')

    # dispensing
    disp = Fluent.dispense()
    disp.RackLabel = df_map.ix[i,'TECAN_dest_labware_name']
    disp.RackType = df_map.ix[i,'TECAN_dest_labware_type']
    disp.Position = df_map.ix[i,'TECAN_dest_target_position']
    disp.Volume = primer_plate_volume
    asp.LiquidClass = liq_cls
    outFH.write(disp.cmd() + '\n')
    
    # tip to waste
    outFH.write('W;\n')
    
    # tracking labware & tip usage
    lw_tracker.add(asp)
    lw_tracker.add(disp, add_tip=False)
    
        
def pip_primers(df_map, outFH, fp_volume=0, rp_volume=0, fp_tube=0, rp_tube=0, 
                liq_cls='Water Free Single', lw_tracker=None):
    """Commands for aliquoting primers
    """
    outFH.write('C;Primers\n')
    primer_plate_volume = 0
    # pipetting non-barcoded primers
    ## forward primer
    if fp_tube > 0 and fp_volume > 0:
        pip_nonbarcode_primer(df_map, outFH, fp_volume, fp_tube, 
                              lw_tracker=lw_tracker)
    else:
        primer_plate_volume += fp_volume
    ## reverse primer
    if rp_tube > 0 and rp_volume > 0:
        pip_nonbarcode_primer(df_map, outFH, rp_volume, rp_tube, 
                              lw_tracker=lw_tracker)
    else:
        primer_plate_volume += rp_volume

    # single or dual barcodes?
    #req_cols_s = ['TECAN_primer_labware', 'TECAN_primer_location']
    #req_cols_d = ['TECAN_primer_labware_F', 'TECAN_primer_location_F',
    #             'TECAN_primer_labware_R', 'TECAN_primer_location_R']
    #barcode_type = None
    #if all(x in df_map.columns.values for x in req_cols_s):
    #    barcode_type = 'single'
    #elif all(x in df_map.columns.values for x in req_cols_d):
    #    barcode_type = 'dual'
    #else:
    #    raise ValueError('Cannot find required TECAN_primer_* columns')
            
    # pipetting barcoded primers
    outFH.write('C;Barcoded primers\n')
    for i in range(df_map.shape[0]):
        pip_primer(i, outFH, df_map,
                   'TECAN_primer_labware_name', 
                   'TECAN_primer_labware_type',
                   'TECAN_primer_target_position',
                   primer_plate_volume, liq_cls, 
                   lw_tracker=lw_tracker)

    ## for each Sample-PCR_rxn_rep, write out asp/dispense commands
    # if barcode_type == 'single':
    #     for i in range(df_map.shape[0]):
    #         pip_primer(i, outFH, df_map,
    #                    'TECAN_primer_labware', 'TECAN_primer_location',
    #                    primer_plate_volume, liq_cls)
    # elif barcode_type == 'dual':
    #     for i in range(df_map.shape[0]):
    #         pip_primer(i, outFH, df_map,
    #                'TECAN_primer_labware_F', 'TECAN_primer_location_F',
    #                primer_plate_volume, liq_cls)
    #     for i in range(df_map.shape[0]):
    #         pip_primer(i, outFH, df_map,
    #                    'TECAN_primer_labware_R', 'TECAN_primer_location_R',
    #                    primer_plate_volume, liq_cls)            
    # else:
    #     raise ValueError('barcode_type not recognized')
                

def pip_samples(df_map, outFH, liq_cls='Water Free Single', lw_tracker=None):
    """Commands for aliquoting samples to each PCR rxn
    """
    outFH.write('C;Samples\n')
    # for each Sample-PCR_rxn_rep, write out asp/dispense commands
    for i in range(df_map.shape[0]):
        # aspiration
        asp = Fluent.aspirate()
        asp.RackLabel = df_map.ix[i,'TECAN_sample_labware_name']
        asp.RackType = df_map.ix[i,'TECAN_sample_labware_type']
        asp.Position = df_map.ix[i,'TECAN_sample_target_position']
        asp.Volume = df_map.ix[i,'TECAN_sample_rxn_volume']
        asp.LiquidClass = liq_cls
        outFH.write(asp.cmd() + '\n')

        # dispensing
        disp = Fluent.dispense()
        disp.RackLabel = df_map.ix[i,'TECAN_dest_labware_name']
        disp.RackType = df_map.ix[i,'TECAN_dest_labware_type']
        disp.Position = df_map.ix[i,'TECAN_dest_target_position']
        disp.Volume = df_map.ix[i,'TECAN_sample_rxn_volume']
        disp.LiquidClass = liq_cls
        outFH.write(disp.cmd() + '\n')

        # tip to waste
        outFH.write('W;\n')

        # tracking labware & tip usage
        lw_tracker.add(asp)
        lw_tracker.add(disp, add_tip=False)


def pip_water(df_map, outFH, pcr_volume=25.0, mm_volume=13.1,
              fp_volume=2.0, rp_volume=2.0, liq_cls='Water Free Single',
              lw_tracker=None):
    """Commands for aliquoting water to each PCR rxn
    """
    outFH.write('C;Water\n')
    # calculate the amount of water
    water_volume = []
    for i in range(df_map.shape[0]):
        samp_volume = df_map.ix[i,'TECAN_sample_rxn_volume']
        w_need = pcr_volume - (samp_volume + mm_volume + fp_volume + rp_volume)
        assert w_need >= 0, 'Water volume is negative: {}'.format(w_need)
        water_volume.append(w_need)
    df_map['TECAN_water_rxn_volume'] = water_volume
    water_target_position = lw_tracker.next_target_position('1.5ml Eppendorf')
    
    # for each Sample-PCR_rxn_rep, write out asp/dispense commands
    for i in range(df_map.shape[0]):
        # aspiration
        asp = Fluent.aspirate()
        #asp.RackLabel = df_map.ix[i,'TECAN_sample_labware_name']
        #asp.RackType = df_map.ix[i,'TECAN_sample_labware_type']
        #asp.Position = df_map.ix[i,'TECAN_sample_target_position']
        asp.RackLabel = 'PCR water tube'
        asp.RackType = '1.5ml Eppendorf'
        asp.Position = water_target_position
        asp.Volume = water_volume[i]
        asp.LiquidClass = liq_cls
        outFH.write(asp.cmd() + '\n')

        # dispensing
        disp = Fluent.dispense()
        disp.RackLabel = df_map.ix[i,'TECAN_dest_labware_name']
        disp.RackType = df_map.ix[i,'TECAN_dest_labware_type']
        disp.Position = df_map.ix[i,'TECAN_dest_target_position']
        disp.LiquidClass = liq_cls
        disp.Volume = water_volume[i]
        outFH.write(disp.cmd() + '\n')

        # tip to waste
        outFH.write('W;\n')

        # tracking labware & tip usage
        lw_tracker.add(asp)
        lw_tracker.add(disp, add_tip=False)

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

def write_report(df_map, outFH, pcr_volume, mm_volume,
                 fp_tube, fp_volume, rp_tube, rp_volume,
                 n_rxn_reps, error_perc=10.0):
    """Writing a report on
    """
    # calculating total volumes
    n_rxn = df_map.shape[0]
    ## total PCR
    total_pcr_volume = pcr_volume * n_rxn
    ## total mastermix
    total_mm_volume = pcr_volume * n_rxn
    ## total primer
    if fp_tube > 0 and fp_volume > 0:
        total_fp_volume = fp_volume * n_rxn
    else:
        total_fp_volume = None
    if rp_tube > 0 and rp_volume > 0:
        total_rp_volume = rp_volume * n_rxn
    else:
        total_rp_volume = None
    ## total water
    total_water_volume = sum(df_map['TECAN_water_rxn_volume'])

    # report
    # number of samples
    outFH.write('# PCR REPORT\n')
    outFH.write('Number of PCR replicates:\t{}\n'.format(n_rxn_reps))
    outFH.write('Number of total PCRs:\t{}\n'.format(n_rxn))
    ## rxn volumes
    outFH.write('# Volumes per PCR (ul)\n')
    write_report_line(outFH, 'Total', pcr_volume)
    write_report_line(outFH, 'Master Mix', mm_volume)
    write_report_line(outFH, 'Forward Primer', fp_volume)
    write_report_line(outFH, 'Reverse Primer', rp_volume)
    ## raw total volumes
    outFH.write('# Total volumes (ul)\n')
    write_report_line(outFH, 'Master Mix', total_mm_volume)
    write_report_line(outFH, 'Forward Primer', total_fp_volume)
    write_report_line(outFH, 'Reverse Primer', total_rp_volume)
    write_report_line(outFH, 'Water', total_water_volume)
    ## total volumes with error
    outFH.write('# Total volumes with (ul; {}% error)\n'.format(error_perc))
    write_report_line(outFH, 'Master Mix', total_mm_volume, error_perc=error_perc)
    write_report_line(outFH, 'Forward Primer', total_fp_volume, error_perc=error_perc)
    write_report_line(outFH, 'Reverse Primer', total_rp_volume, error_perc=error_perc)
    write_report_line(outFH, 'Water', total_water_volume, error_perc=error_perc)
    # samples
    outFH.write('')


# main
if __name__ == '__main__':
    pass


