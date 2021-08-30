#!/usr/bin/env python
from __future__ import division, print_function
from astropy.io import fits
from astropy.wcs import WCS
import csv
from astropy.io.votable import parse
from argparse import ArgumentParser
import glob

def main(args):
    


if __name__ == "__main__":
    parser = ArgumentParser("makeVOTable", description="combines all the detections and saves all the data in a vo table")
    parser.add_argument("--obs",required=True, type=int, help="the obs id")
    args = parser.parse_args()

    main(args)
