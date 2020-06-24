from __future__ import division
from __future__ import print_function
from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
from astropy.nddata Cutout2D
import os.path
from argparse import ArgumentParser
import time
from datetime import datetime, timedelta
import csv
from subprocess import call


def main(obs, freq, timeStep, prefix):
    

if __name__ == "__main__":
    parser = ArguementParser("ORB_INFO", description="obtains all freq and position info for the pass")
    parser.add_argument("--obs", required=True, help="the obs id")
    parser.add_argument("--freq",default=768, type=int, help="the number of frequency channels to process")
    parser.add_arguement("--timeStep", required=True, tyep=int, help="the time step at which sat info is being searched for")
    parser.add_argument("--prefix", required=True,help="the prefix used in the output files")
    parser.add_argument("--debug",type=bool,defaut=False, help="runs program in debug mode")
    args = parser.parser_args()
    global debug
    debug = args.debug

    if debug:
        print("running Orb_info in debug mode")

    main(args.obs, args.freq, args.timeStep, args.prefix)


