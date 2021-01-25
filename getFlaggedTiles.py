#!/usr/bin/env python
from __future__ import print_function, division
from argparse import ArgumentParser
from astropy.io import fits
import numpy as np


def main(args):
    hdu = fits.open(args.metafits)
    flags = hdu[1].data["Flag"].astype(bool)
    flaggedTiles = set(hdu[1].data["Antenna"][flags == True])

    for tile in flaggedTiles:
        print(tile)


if __name__ == "__main__":
    parser = ArgumentParser("getFlaggedTiles", description="obtains the flagged tiles from the .metafits")
    parser.add_argument("--metafits", required=True, help="the metafits file")
    args = parser.parse_args()

    main(args)