#!/usr/bin/evn python
from __future__ import print_function, division
from argparse import ArgumentParser
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

def main(args):

    hdu = fits.open(args.metafits)
    pointing = SkyCoord(hdu[0].header['RA'], hdu[0].header['DEC'], unit=(u.degree, u.degree))

    print(pointing.to_string('hmsdms'))

if __name__ == "__main__":
    parser = ArgumentParser('pointing', description="gets the pointing direction of the obs")
    parser.add_argument("--metafits", required=True, help="the metafits file")
    args = parser.parse_args()

    main(args)
