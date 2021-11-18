from __future__ import print_function
# import
## batteries
import os
import re
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
    desc = 'Convert a mapping file to a NGS amplicon worklist file for the TECAN robot'
    return desc

def parse_args(test_args=None, subparsers=None):
    desc = get_desc()
    epi = """DESCRIPTION:
    Convert a QIIME-formatted mapping file to a GWL file, which is used by the TECAN
    robot to conduct the NGS amplicon PCR prep (ie., combining MasterMix, primers, samples, etc).
    The extra columns in the mapping file designate the SOURCE of samples and primers;
    the DESTINATION (plate & well) is set by this script. 

    EXTRA COLUMNS in MAPPING FILE:

    * "TECAN_sample_labware_name" = The sample labware name on the robot worktable. Whatever name you want to use! 
    * "TECAN_sample_labware_type" = The type of labware containing samples (eg., '96 Well Eppendorf TwinTec PCR')
    * "TECAN_sample_target_position" = The well or tube location (a number)
    * "TECAN_sample_rxn_volume" = The volume of sample to use per PCR (ul)
    * "TECAN_primer_labware_name" = The primer plate labware name on the robot worktable
    * "TECAN_primer_labware_type" = The primer plate labware type on the robot worktable
    * "TECAN_primer_target_position" = The well location (1-96 or 1-384)

    PRIMERS:
    * For single-indexed primers (eg., EMP primers), the non-barcoded primer should be pre-added to either
      the mastermix (adjust the volume used for this script!) or each of the barcoded primers.
    * If --prm-volume set to 0 or --prm-in-mm, then primers are skipped. 

    CONTROLS:
    * For the positive & negative controls, include them in the mapping file.
    * If the controls (or samples) are provided in a tube, use "1.5ml Eppendorf" 
      for the "TECAN_sample_labware_type" column. The labware name can be whatever you want.

    LABWARE:
    * In order to use a plate adapter, use "PCR Adapter 96 Well and 96 Well Eppendorf TwinTec PCR"

    WATER:
    * If  mastermix + primer + sample = total_reaction_volume, then no water is aliquoted
    * If --water-in-mm, then water is skipped

    MISC NOTES:
    * By default, water is assumed to be in the mastermix, and a 25ml trough will be used to hold the mastermix
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
    dest.add_argument('--dest-name', type=str, default='Destination plate',
                      help='Distination labware name (default: %(default)s)')
    dest.add_argument('--dest-type', type=str, default='96 Well Eppendorf TwinTec PCR',
                      choices=['96 Well Eppendorf TwinTec PCR',
                               'PCR Adapter 96 Well and 96 Well Eppendorf TwinTec PCR',
                               '384 Well Biorad PCR',
                               'PCR Adapter 384 Well and 384 Well Biorad PCR'],
                      help='Destination labware type (default: %(default)s)')
    dest.add_argument('--dest-start', type=int, default=1,
                      help='Start well number on destination plate (default: %(default)s)')
    dest.add_argument('--rxns', type=int, default=3,
                      help='Number of replicate PCRs per sample (default: %(default)s)')
    
    ## Reagents
    rgnt = parser.add_argument_group('Reagents')
    rgnt.add_argument('--pcr-volume', type=float, default=25.0,
                        help='Total volume per PCR (default: %(default)s)')
    rgnt.add_argument('--mm-volume', type=float, default=13.1,
                      help='MasterMix volume per PCR (default: %(default)s)')
    rgnt.add_argument('--prm-volume', type=float, default=2.0,
                         help='Primer volume per PCR (default: %(default)s)')    
    rgnt.add_argument('--error-perc', type=float, default=10.0,
                        help='Percent of extra total reagent volume to include (default: %(default)s)')
    rgnt.add_argument('--mm-labware-type', type=str, default='25ml_1 waste',
                      help='Labware type for mastermix (default: %(default)s)')
    rgnt.add_argument('--mm-one-source', action='store_true', default=False,
                      help='Just one mastermix labware instead of one per destination plate? (default: %(default)s)')
    rgnt.add_argument('--prm-in-mm', action='store_true', default=False,
                         help='Primer pre-added to mastermix? (default: %(default)s)')
    rgnt.add_argument('--water-in-mm', action='store_true', default=False,
                         help='Water pre-added to mastermix? (default: %(default)s)')
    
    # Liquid classes
    liq = parser.add_argument_group('Liquid classes')
    liq.add_argument('--mm-liq', type=str, default='MasterMix Free Single Wall Disp',
                      help='Mastermix liquid class (default: %(default)s)')
    liq.add_argument('--primer-liq', type=str, default='Water Free Single Wall Disp',
                      help='Primer liquid class (default: %(default)s)')
    liq.add_argument('--sample-liq', type=str, default='Water Free Single Wall Disp',
                      help='Sample liquid class (default: %(default)s)')
    liq.add_argument('--water-liq', type=str, default='Water Free Single Wall Disp',
                      help='Water liquid class (default: %(default)s)')
    liq.add_argument('--n-tip-reuse', type=int, default=4,
                     help='Number of tip reuses for applicable reagents (default: %(default)s)')
    liq.add_argument('--n-multi-disp', type=int, default=1,
                     help='Number of tip reuses for applicable reagents (default: %(default)s)')
    
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
    df_map = add_dest(df_map,
                      dest_labware=args.dest_name,
                      dest_type=args.dest_type,
                      dest_start=args.dest_start,
                      rxn_reps=args.rxns)
    df_map = check_rack_labels(df_map)
    
    # gwl construction
    TipTypes = ['FCA, 1000ul SBS', 'FCA, 200ul SBS',
                'FCA, 50ul SBS', 'FCA, 10ul SBS']     
    gwl = Fluent.gwl(TipTypes)

    # Reordering dest if plate type is 384-well
    df_map = Utils.reorder_384well(df_map, gwl,
                                   labware_name_col='TECAN_dest_labware_name',
                                   labware_type_col='TECAN_dest_labware_type',
                                   position_col='TECAN_dest_target_position')

        
    ## mastermix
    pip_mastermix(df_map, gwl,
                  mm_labware_type=args.mm_labware_type,
                  mm_volume=args.mm_volume, 
                  liq_cls=args.mm_liq,
                  n_tip_reuse=args.n_tip_reuse,
                  n_multi_disp=args.n_multi_disp,
                  mm_one_source=args.mm_one_source)

    ## primers
    if args.prm_in_mm == True:
        args.prm_volume = 0
    if args.prm_volume > 0:
        pip_primers(df_map, gwl,
                    prm_volume=args.prm_volume,
                    liq_cls=args.primer_liq)
    else:
        msg = 'WARNING: primers skipped; make sure that primers are added to the mastermix!'
        print(msg, file=sys.stderr)

    ## samples
    pip_samples(df_map, gwl, liq_cls=args.sample_liq)
    
    ## water
    df_map = calc_water_needed(df_map,
                               pcr_volume=args.pcr_volume,
                               mm_volume=args.mm_volume,
                               prm_volume=args.prm_volume,
                               water_in_mm=args.water_in_mm)
    if sum(df_map['TECAN_water_rxn_volume']) > 0:
        pip_water(df_map, gwl, liq_cls=args.water_liq)
    else:
        msg = 'WARNING: water skipped; make sure that water is added to the mastermix!'
        print(msg, file=sys.stderr)
    
    ## writing out worklist (gwl) file
    gwl_file = args.prefix + '.gwl'
    gwl.write(gwl_file)
    
    # Report (total volumes; sample truncation; samples)
    report_file = args.prefix + '_report.txt'
    with open(report_file, 'w') as repFH:
        write_report(df_map, outFH=repFH,
                     pcr_volume=args.pcr_volume,
                     mm_volume=args.mm_volume,
                     prm_volume=args.prm_volume,
                     n_rxn_reps=args.rxns,
                     error_perc=args.error_perc)
        
    # making labware table
    lw = Labware.labware()
    lw.add_gwl(gwl)
    lw_df = lw.table()
    lw_file = args.prefix + '_labware.txt'
    lw_df.to_csv(lw_file, sep='\t', index=False)
        
    # Mapping file with destinations
    df_file = args.prefix + '_map.txt'
    df_map['TECAN_water_rxn_volume'] = df_map['TECAN_water_rxn_volume'].round(2)
    df_map['TECAN_dest_target_position'] = df_map['TECAN_dest_target_position'].astype(int)
    df_map['TECAN_pcr_rxn_rep'] = df_map['TECAN_pcr_rxn_rep'].astype(int)
    df_map.to_csv(df_file, sep='\t', index=False, na_rep='NA')

    # Plate map file for Bio-Rad PrimePCR software (designates: sampleID <--> wellID)
    biorad_files = PrimerPCR_plate_map(df_map, prefix=args.prefix)
    
    # status on files written
    Utils.file_written(gwl_file)
    Utils.file_written(lw_file)
    Utils.file_written(report_file)
    [Utils.file_written(x) for x in biorad_files]
    Utils.file_written(df_file)
    
    # Return
    return (gwl_file, report_file, df_file, lw_file)
        
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
    ## mastermix
    if args.mm_volume > args.pcr_volume:
        msg = 'MasterMix volume > total PCR volume'
        raise ValueError(msg)
    if args.mm_volume / 2.0 > args.pcr_volume:
        msg = 'WARNING: MasterMix volume > half of PCR volume'
        print(msg, file=sys.stderr)
    ## primers
    assert args.prm_volume >= 0.0, 'Primer volume must be >= 0'
    assert args.pcr_volume > 0.0, 'PCR volume must be > 0'    
        
def map2df(mapfile, row_select=None):
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
        msg = 'WARNING: Duplicated sample values in the mapping file. Did intend for this?'
        print(msg, file=sys.stderr)

    # checking for unique barcode locations (NaN's filtered out)
    dups = df_map[req_cols].dropna().duplicated(keep=False)
    if any(dups):
        msg = 'WARNING: Duplicated barcodes in the mapping file!'
        print(msg, file=sys.stderr)
        
    # checking sample volumes
    msg = 'WARNING: sample volume > mastermix volume'
    for sv in df_map['TECAN_sample_rxn_volume']:
        if sv > args.mm_volume:
            print(msg, file=sys.stderr)
        if sv < 0:
            raise ValueError('Sample volume < 0')

    # sampleID rename
    if '#SampleID' in df_map.columns.values:
        df_map.rename(columns={'#SampleID':'SampleID'}, inplace=True)
        
    # making sure labware names are "TECAN worklist friendly"
    df_map = Utils.rm_special_chars(df_map, 'TECAN_sample_labware_name')
    df_map = Utils.rm_special_chars(df_map, 'TECAN_primer_labware_name')

    # removing "tube" from end of labware type (if present)
    Utils.rm_tube(df_map, 'TECAN_sample_labware_type')
    Utils.rm_tube(df_map, 'TECAN_primer_labware_type')
        
def check_rack_labels(df_map):
    """Removing '.' for rack labels (causes execution failures)
    """
    cols = ['TECAN_sample_labware_name', 'TECAN_primer_labware_name', 'TECAN_dest_labware_name']
    for x in cols:
        df_map[x] = [y.replace('.', '_') for y in df_map[x].tolist()]
    return df_map
        
def add_dest(df_map, dest_labware,
             dest_type='96 Well Eppendorf TwinTec PCR',
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
    db = Fluent.db()    
    dest_start= int(dest_start)
    
    # labware type found in DB?
    db.get_labware(dest_type)
    positions = db.get_labware_wells(dest_type)

    # init destination df
    sample_col = df_map.columns[0]
    cols = [sample_col, 'TECAN_pcr_rxn_rep',
            'TECAN_dest_labware_name',
            'TECAN_dest_labware_type', 
            'TECAN_dest_target_position']    
    ncol = len(cols)
    nrow = df_map.shape[0] * rxn_reps        # number of rxns
    df_dest = pd.DataFrame(np.nan, index=range(nrow), columns=cols)

    # number of destination plates required
    n_dest_plates = round(df_dest.shape[0] / positions + 0.5, 0)
    if n_dest_plates > 1:
        msg = ('WARNING: Not enough wells for the number of samples.' 
        ' Using multiple destination plates')
        print(msg, file=sys.stderr)
        
    
    # filling destination df
    orig_dest_labware = dest_labware
    for i,(sample,rep) in enumerate(product(df_map.iloc[:,0], range(rxn_reps))):
        # dest location
        dest_position = i + dest_start
        dest_position = positions if dest_position % positions == 0 else dest_position % positions 

        # destination plate name
        if n_dest_plates > 1:
            x = round((i + dest_start) / positions + 0.499, 0)        
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

def pip_mastermix(df_map, gwl, mm_labware_type='25ml_1 waste',
                  mm_volume=13.1, n_tip_reuse=6, n_multi_disp=4,
                  liq_cls='MasterMix Free Multi', mm_one_source=False):
    """Writing worklist commands for aliquoting mastermix.
    Using 1-asp-multi-disp with 200 ul tips.
    Method:
    * calc max multi-dispense for 200 ul tips & mm_volume
    * for 1:n_dispense
      * determine how many disp left for channel (every 8th)
      * if n_disp_left < n_disp: n_disp = n_disp_left
      * calc total volume: n_disp * mm_volume
    """
    gwl.add(Fluent.Comment('MasterMix'))

    # copying df
    df = df_map.copy()
    
    # separate regent dispense per designation plate
    cols = ['TECAN_dest_labware_name', 'TECAN_dest_labware_type']
    df_f = df.loc[:,cols].drop_duplicates()
    df_f.reset_index(inplace=True)
    func = lambda row: gwl.db.get_labware_wells(row['TECAN_dest_labware_type'])
    df_f['wells'] = df_f.apply(func, axis=1)

    ## ordering df for proper tip reuse
    if n_multi_disp == 1:
        x = cycle(range(8))
        df['CHANNEL_ORDER'] = [next(x) for y in range(df.shape[0])]
        x = cycle(range(n_tip_reuse))
        df['TIP_BATCH'] = Utils.tip_batch(df['CHANNEL_ORDER'], n_tip_reuse)
        df.sort_values(by=['TIP_BATCH',
                           'CHANNEL_ORDER',
                           'TECAN_dest_target_position'], inplace=True)
        df.reset_index(inplace=True)
        
    # dispense
    for i in range(df_f.shape[0]):
        # all records for 1 plate
        RackLabel = df_f.loc[i,'TECAN_dest_labware_name']
        df_tmp = df.loc[df['TECAN_dest_labware_name'] == RackLabel,]
        df_tmp.reset_index(inplace=True)
        # dispense single or with reagent distribution
        if n_multi_disp == 1:
            # creating asp-dispense
            liq_cls = re.sub('Multi', 'Single', liq_cls)
            for ii in range(df_tmp.shape[0]):
                # aspiration
                asp = Fluent.Aspirate()
                if mm_one_source == True:
                    asp.RackLabel = 'Mastermix'
                else:
                    asp.RackLabel = 'Mastermix[{0:0>3}]'.format(i + 1)
                asp.RackType = mm_labware_type
                asp.Position = 1
                asp.Volume = mm_volume
                asp.LiquidClass = liq_cls
                gwl.add(asp)

                # dispensing
                disp = Fluent.Dispense()
                disp.RackLabel = df_tmp.loc[ii,'TECAN_dest_labware_name']
                disp.RackType = df_tmp.loc[ii,'TECAN_dest_labware_type']
                disp.Position = df_tmp.loc[ii,'TECAN_dest_target_position']
                disp.Volume = mm_volume
                disp.LiquidClass = liq_cls
                gwl.add(disp)
                
                # waste
                if (ii + 1) % n_tip_reuse == 0 or ii + 1 == df_tmp.shape[0]:
                    gwl.add(Fluent.Waste())
                
        # using reagent distribution
        else:
            # finding positions to exclude
            all_wells = [x+1 for x in range(df_f.loc[i,'wells'])]
            target_pos = df_tmp['TECAN_dest_target_position'].tolist()
            to_exclude = set(all_wells) - set(target_pos)            
            # creating reagnet distribution command
            rd = Fluent.Reagent_distribution()
            if mm_one_source == True:
                rd.SrcRackLabel = 'Mastermix'
            else:
                rd.SrcRackLabel = 'Mastermix[{0:0>3}]'.format(i + 1)
            rd.SrcRackType = mm_labware_type
            rd.SrcPosStart = 1
            rd.SrcPosEnd = 1
            # dispense parameters
            rd.DestRackLabel = df_f.loc[i,'TECAN_dest_labware_name']
            rd.DestRackType = df_f.loc[i,'TECAN_dest_labware_type']
            rd.DestPosStart = 1
            rd.DestPosEnd = df_f.loc[i,'wells']
            # other
            rd.Volume = mm_volume
            rd.LiquidClass = liq_cls
            rd.NoOfDiTiReuses = n_tip_reuse
            rd.NoOfMultiDisp = n_multi_disp
            rd.Direction = 0
            rd.ExcludedDestWell = ';'.join([str(x) for x in list(to_exclude)])
            # adding to gwl object
            gwl.add(rd)

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
        pip_primer(i, gwl, df_map,
                   'TECAN_primer_labware_name', 
                   'TECAN_primer_labware_type',
                   'TECAN_primer_target_position',
                   prm_volume, liq_cls)
        
    # adding break
    gwl.add(Fluent.Break())

                
def pip_samples(df_map, gwl, liq_cls='Water Free Single'):
    """Commands for aliquoting samples to each PCR rxn
    """
    gwl.add(Fluent.Comment('Samples'))
    # for each Sample-PCR_rxn_rep, write out asp/dispense commands
    for i in range(df_map.shape[0]):
        # aspiration
        asp = Fluent.Aspirate()
        asp.RackLabel = df_map.loc[i,'TECAN_sample_labware_name']
        asp.RackType = df_map.loc[i,'TECAN_sample_labware_type']
        asp.Position = df_map.loc[i,'TECAN_sample_target_position']
        asp.Volume = round(df_map.loc[i,'TECAN_sample_rxn_volume'], 1)
        asp.LiquidClass = liq_cls
        gwl.add(asp)

        # dispensing
        disp = Fluent.Dispense()
        disp.RackLabel = df_map.loc[i,'TECAN_dest_labware_name']
        disp.RackType = df_map.loc[i,'TECAN_dest_labware_type']
        disp.Position = df_map.loc[i,'TECAN_dest_target_position']
        disp.Volume = round(df_map.loc[i,'TECAN_sample_rxn_volume'], 1)
        disp.LiquidClass = liq_cls
        gwl.add(disp)

        # waste
        gwl.add(Fluent.Waste())
        
    # adding break
    gwl.add(Fluent.Break())

def calc_water_needed(df_map, pcr_volume=25.0, mm_volume=13.1,
                      prm_volume=2.0, water_in_mm=False):
    """Calculating the amount of water to reach full rxn volume
    """
    if water_in_mm == True:
        df_map['TECAN_water_rxn_volume'] = 0
    else:
        water_volume = []
        for i in range(df_map.shape[0]):
            samp_volume = df_map.loc[i,'TECAN_sample_rxn_volume']
            w_need = pcr_volume - (samp_volume + mm_volume + prm_volume)
            if w_need <= 0:
                w_need = 0
                #assert w_need >= 0, 'Water volume is negative: {}'.format(w_need)
            water_volume.append(w_need)
        df_map['TECAN_water_rxn_volume'] = water_volume
    return df_map
        
def pip_water(df_map, gwl, liq_cls='Water Free Single'):
    """Commands for aliquoting water to each PCR rxn
    """
    gwl.add(Fluent.Comment('Water'))
    
    # for each Sample-PCR_rxn_rep, write out asp/dispense commands
    for i in range(df_map.shape[0]):
        if df_map.loc[i,'TECAN_water_rxn_volume'] <= 0:
            continue
        # aspiration
        asp = Fluent.Aspirate()
        asp.RackLabel = '25ml_1[001]'
        asp.RackType = '25ml_1 waste'
        asp.Position = 1 
        asp.Volume = round(df_map.loc[i,'TECAN_water_rxn_volume'], 1)
        asp.LiquidClass = liq_cls
        gwl.add(asp)

        # dispensing
        disp = Fluent.Dispense()
        disp.RackLabel = df_map.loc[i,'TECAN_dest_labware_name']
        disp.RackType = df_map.loc[i,'TECAN_dest_labware_type']
        disp.Position = df_map.loc[i,'TECAN_dest_target_position']
        disp.LiquidClass = liq_cls
        disp.Volume = round(df_map.loc[i,'TECAN_water_rxn_volume'], 1)
        gwl.add(disp)

        ## waste
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
            biorad_file = prefix + '_BIORAD-{}.csv'.format(labware.replace(' ', '_'))
            df_biorad.loc[labware].to_csv(biorad_file, sep=sep, index=False, na_rep='')
            biorad_files.append(biorad_file)
    else:
        biorad_file = prefix + '_BIORAD.csv'
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

def write_report(df_map, outFH, pcr_volume, mm_volume,
                 prm_volume, n_rxn_reps, error_perc=10.0):
    """Writing a report on
    """    
    # calculating total volumes
    n_rxn = df_map.shape[0]
    ## total PCR
    total_pcr_volume = pcr_volume * n_rxn
    ## total mastermix
    total_mm_volume = mm_volume * n_rxn
    ## total primers
    total_prm_volume = prm_volume * n_rxn
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
    write_report_line(outFH, 'Primers', prm_volume)
    ## raw total volumes
    outFH.write('# Total volumes (ul)\n')
    write_report_line(outFH, 'Master Mix', total_mm_volume)
    write_report_line(outFH, 'Water', total_water_volume)
    write_report_line(outFH, 'Primers', total_prm_volume)
    ## total volumes with error
    outFH.write('# Total volumes with (ul; {}% error)\n'.format(error_perc))
    write_report_line(outFH, 'Master Mix', total_mm_volume, error_perc=error_perc)
    write_report_line(outFH, 'Water', total_water_volume, error_perc=error_perc)
    write_report_line(outFH, 'Primers', total_prm_volume, error_perc=error_perc)
    # samples
    outFH.write('')


# main
if __name__ == '__main__':
    pass


