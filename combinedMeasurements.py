#!/usr/bin/env python
from __future__ import division, print_function
import csv
from astropy.io import fits
from astropy.wcs import WCS
from argparse import ArgumentParser
import csv
from astropy.table import QTable
import astropy.units as u
import numpy as np
from tqdm import tqdm

def main(args):
    
    timeSteps = range(args.t1, args.t2 + 1)
    x_array, y_array, ra_array, dec_array, az_array, el_array\
     = [], [], [], [], [], []

    snr_array, source_array, flux_array = [], [], []

    utc_array, freq_array, timeStep_array, channel_array\
    = [], [], [], []
    head_array, tail_array, floodfill_array, seed_array=\
    [], [], [], []

    for t in tqdm(timeSteps):

        with open("measuredValues_t" + str(t).zfill(4) + ".csv") as csv_file:

            csv_reader = csv.reader(csv_file, delimiter=",")

            for row in csv_reader:

                ra, dec, snr, sourceNo, flux, az, el, freq,\
                utc, channel, ts, x, y, head, tail, flood, seed = row

                if str(x) == "x":
                    continue

                else:
                    x_array.append(float(x))
                    y_array.append(float(y))
                    ra_array.append(float(ra))
                    dec_array.append(float(dec))
                    az_array.append(float(az))
                    el_array.append(float(el))
                    snr_array.append(float(snr))
                    source_array.append(float(sourceNo))
                    flux_array.append(float(flux))
                    utc_array.append(utc)
                    freq_array.append(float(freq))
                    timeStep_array.append(float(ts))
                    channel_array.append(float(channel))
                    head_array.append(head)
                    tail_array.append(tail)
                    floodfill_array.append(flood)
                    seed_array.append(seed)


    hdu = fits.open(str(args.prefix) + "Sigma" + "RFIBinaryMapPeakFlux-t0001.fits")
    header = hdu[0].header


    ra_array = ra_array * u.deg
    dec_array = dec_array * u.deg
    timeStep_array = timeStep_array * u.s
    tab = QTable([ra_array, dec_array, snr_array, source_array, flux_array, az_array, el_array,\
                 freq_array, utc_array, channel_array, timeStep_array, x_array, y_array,
                 head_array, tail_array, floodfill_array, seed_array],
                 names=("ra","dec","SNR","sourceNo","FluxDensity","az","elv","freq","utc"\
                 ,"channelNo","timeStep","x","y","headFile","tailFile","floodfillSigma","seedSigma"),
                 meta={"CTYPE1" : header["CTYPE1"],
                       "CTYPE2" : header["CTYPE2"],
                       "CTYPE3" : header["CTYPE3"],
                       "CTYPE4" : header["CTYPE4"],
                       "CRVAL1" : header["CRVAL1"],
                       "CRVAL2" : header["CRVAL2"],
                       "CRVAL3" : header["CRVAL3"],
                       "CRVAL4" : header["CRVAL4"],
                       "CRPIX1" : header["CRPIX1"],
                       "CRPIX2" : header["CRPIX2"],
                       "CRPIX3" : header["CRPIX3"],
                       "CRPIX4" : header["CRPIX4"],
                       "CDELT1" : header["CDELT1"],
                       "CDELT2" : header["CDELT2"],
                       "CDELT3" : header["CDELT3"],
                       "CDELT4" : header["CDELT4"],
                       "NAXIS" : header["NAXIS"],
                       "NAXIS1" : header["NAXIS1"],
                       "NAXIS2" : header["NAXIS2"]})

    tab.write(str(args.obs) + "-"+ str(args.hpc)+"-measurements.fits")




    

if __name__ == "__main__":
    parser = ArgumentParser("combinedMeasurements", description="combines all the csv files")
    parser.add_argument("--t1", default=1, type=int, help="the start timeStep")
    parser.add_argument("--t2", default=55, type=int, help="the last timeStep")
    parser.add_argument("--obs", required=True, type=int, help="the obs id")
    parser.add_argument("--hpc",default="pawsey", help="the name of the hpc used to process data")
    parser.add_argument("--prefix", required=True, help="the file prefix to obtain the wcs info")
    args = parser.parse_args()

    main(args)
    
