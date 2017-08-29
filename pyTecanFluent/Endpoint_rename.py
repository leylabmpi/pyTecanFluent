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


# functions
def get_desc():
    desc = 'Renaming an BioRad Endpoint file based on PCR mapping file'
    return desc

def parse_args(test_args=None, subparsers=None):
    # desc
    desc = get_desc()
    epi = """DESCRIPTION:
    Renaming an BioRad Endpoint file based on PCR mapping file.
    
    The input is 2 files: 
    * Endpoint file (export all from BioRad software)
    * PCR mapping file which has the names & locations of each sample in the PCR plate. 
      * Names column: "#SampleID" 
      * Location (well) column: "TECAN_dest_location"

    #* Sample labware  (eg., "96 Well[001]")
    #* Sample location (numeric value; minimum of 1)
    #* Sample concentration (numeric value; units=ng/ul)
    
    #Notes:
    #* You can designate the input table columns for each value (see options).
    #* Sample locations in plates numbered are column-wise. 
    #* All volumes are in ul.
    """
    if subparsers:
        parser = subparsers.add_parser('PCR_pool', description=desc, epilog=epi,
                                       formatter_class=argparse.RawTextHelpFormatter)
    else:
        parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                         formatter_class=argparse.RawTextHelpFormatter)

    # args
    ## I/O
    groupIO = parser.add_argument_group('I/O')
    groupIO.add_argument('infile', metavar='infile', type=str,
                         help='An excel or tab-delim file of ')
    groupIO.add_argument('--prefix', type=str, default='TECAN_dilute',
                         help='Output file name prefix (default: %(default)s)')


    ## PCR success/failure
    
    
    ## destination plate
    dest = parser.add_argument_group('Destination plate')
    dest.add_argument('--dest', type=str, default='96 Well[001]',
                      help='Destination plate labware ID on TECAN worktable (default: %(default)s)')
    dest.add_argument('--desttype', type=str, default='96',
                      choices=['96','384'],
                      help='Destination plate labware type (default: %(default)s)')
    dest.add_argument('--deststart', type=int, default=1,
                      help='Start well number on destination OD plate (default: %(default)s)')
#    dest.add_argument('--destlabware', type=str, default='384 Well[004]',
#                      help='Destination Labware ID on the TECAN worktable (default: %(default)s)')
#                      default='96-well:96 Well[008],384-well:384 Well[004]',
#                      help='Choices for the destination labware name base on --desttype')

    # parse & return
    if test_args:
        args = parser.parse_args(test_args)
        return args

    return parser


def main(args=None):
    """
    algorithm
    ---------
    * load 
      * Endpoint file
      * Mapping file from map2robot
    * join tables by location
      * For endpoint file: convert "Well" to location
      * Just 
    * call PCR success/failure/pass based on RFU values
      * Place in "Call" column
      * User defined: "pass" (failures to be sequenced) 
      * Optional: user-defined End RFU cutoff 
    * output
      * write an updated Endpoint file w/ samples & successes/failures/passes
         * also, locations are edited
      * write a truncated mapping file of PCR failures
        * this can be used for subsequent PCR attempts

    """
    # Input
    if args is None:
        args = parse_args()
    check_args(args)


# main
if __name__ == '__main__':
    pass


