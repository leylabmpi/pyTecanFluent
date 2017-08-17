# -*- coding: utf-8 -*-
# import
## batteries
from __future__ import print_function
import os
import sys
import argparse
## 3rd party
from pyTecanFluent import Map2Robot
from pyTecanFluent import Dilute
from pyTecanFluent import QPCR

# main
def main(args=None):
  if args is None:
    args = sys.argv[1:]

  # main parser
  desc = 'pyTecanFluent: a python interface to the TECAN Fluent'
  epi = """DESCRIPTION:
  Create worklist files for the TECAN FluentControl software (*.gwl files)
  """

  parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                   formatter_class=argparse.RawTextHelpFormatter)

  # subparsers
  subparsers = parser.add_subparsers()
  ## map2robot
  map2robot = Map2Robot.parse_args(subparsers=subparsers)
  map2robot.set_defaults(func=Map2Robot.main)
  ## dilute
  dilute = Dilute.parse_args(subparsers=subparsers)
  dilute.set_defaults(func=Dilute.main)
  ## qpcr
  qpcr = QPCR.parse_args(subparsers=subparsers)
  qpcr.set_defaults(func=QPCR.main)

  # parsing args
  if args:
    args = parser.parse_args(args)
  else:
    args = parser.parse_args()

  # 
  if len(vars(args)) > 0:
    args.func(args)
  else:
    parser.parse_args(['--help'])
    
    
if __name__ == '__main__':
    main()
