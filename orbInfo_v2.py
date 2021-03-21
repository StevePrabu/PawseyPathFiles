#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from datetime import datetime, timedelta
from argparse import ArgumentParser
import spacetracktool as st
from spacetracktool import operations as ops
from skyfield.api import EarthSatellite
from skyfield.api import Topos, load
from tqdm import tqdm
from astropy.nddata import Cutout2D
import scipy.optimize as opt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.coordinates import EarthLocation
from astropy.time import Time
from astropy import units as u
import csv
import json
from os import path
from warnings import filterwarnings
filterwarnings("ignore")
import sys
from uncertainties import ufloat
from uncertainties.umath import * 

obs_loc = EarthLocation(lon=116.67083333*u.deg, lat=-26.70331941*u.deg, height=377.827*u.m)

def eprint(*args, **kwargs):
    """
    prints function in stderr so that the print is logged
    """
    print(*args, file=sys.stderr, **kwargs)


def obtainTimeSteps(noradid):
    """
    returns the timeSteps to search for signal
    by searching the json file made by satSearch.py

    Parameters
    ----------
    noradid         : the norad id

    Returns
    -------
    headTimeSteps   : the list of head timeSteps
    tailTimeSteps   : the list of tail timeSteps
    midTimeSteps    : the list of mid timeSteps  

    """
    with open("filtered_summary.json") as f:
        filtered_summary = json.load(f)

    for line in filtered_summary:
        if filtered_summary[line]["norad"] == str(noradid):
            timeSteps = list(map(int,filtered_summary[line]["timeSteps"].split("-")))

    headTimeSteps, tailTimeSteps, midTimeSteps = [], [] , []

    for t in timeSteps:
        headTimeSteps.append(t)
        tailTimeSteps.append(t+1)
    
    midTimeSteps = [i for i in headTimeSteps if i in tailTimeSteps]

    return headTimeSteps, tailTimeSteps, midTimeSteps




def obtainTLE(noradid):
    """
    for a given norad id, obtains the tle for the relevant date range

    Parameters
    ----------
    noradid : the noradid of the satellite 

    Returns
    ------
    returns the tle of the satellite (tle = line1 + line2 + line3)
    """

    time2 = refUTC + timedelta(hours=24)
    day1, month1, year1 = str(refUTC.day).zfill(2), str(refUTC.month).zfill(2), str(refUTC.year)
    day2, month2, year2 = str(time2.day).zfill(2), str(time2.month).zfill(2), str(time2.year)
    date_range = ops.make_range_string(year1 + "-" + month1 + "-" + day1, year2 + "-" + month2 + "-" + day2)

    if path.exists(str(noradid) + ".txt") == False:

        if debug:
            eprint("requesting file from server")

        result = query.tle_query(epoch=date_range,norad_cat_id=noradid)

        with open(str(noradid) + ".txt", "w") as outfile:
            json.dump(result.json(), outfile)

        line1 = result.json()[0]["OBJECT_NAME"]
        line2 = result.json()[0]["TLE_LINE1"]
        line3 = result.json()[0]["TLE_LINE2"]

    else:

        if debug:
            eprint("tle file found on disk. Not downloading.")

        with open(str(noradid) + ".txt") as json_file:
            result = json.load(json_file)

        line1 = result[0]["OBJECT_NAME"]
        line2 = result[0]["TLE_LINE1"]
        line3 = result[0]["TLE_LINE2"]


    return line1 , line2, line3




def getSatMWA(line1, line2, line3):
    """
    creates a satellite object and observer object
    
    Paramters
    ---------
    line1   : tle line 1
    line2   : tle line 2
    line3   : tle line 3

    Returns
    -------
    satellite   : the sat object
    mwa         : the observer object
    ts          : the time object
    """

    ts = load.timescale(builtin=True)
    satellite = EarthSatellite(line2, line3, line1, ts)
    mwa = Topos("26.701276778 S", "116.670846137 E", elevation_m= 377.827)

    return satellite, mwa, ts



def maskData(data,  UT, mwa, sat, wcs, ts, plot=False):
    """
    this functions masks everying outside a 10 km cone from satellite

    Parameters
    ----------
    data    : the output data from RFISeeker
    UT      : the UTC of the timeStep
    mwa     : the observer object
    sat     : the satellite object
    wcs     : world coordinate system
    ts      : timescale object (required by skyfield)

    Returns
    -------
    masked version of the input data variable
    """

    ## propogate the satellite to the utc
    time = ts.utc(UT.year, UT.month,  UT.day, UT.hour, UT.minute, UT.second)
    sat.at(time)
    difference = sat - mwa
    topocentric = difference.at(time)
    ra, dec, distance = topocentric.radec()
    ra = np.degrees(ra.radians)
    dec = np.degrees(dec.radians)
    px, py = wcs.all_world2pix([[ra, dec]], 1)[0]
    LOS_range = distance.m
    
    ## mask the data
    radius = np.degrees(50000/LOS_range) # cone radius (s = r * theta)
    number_of_pixels = radius/imgScale ### check this line

    if plot:
        print("sat location")
        print("x {0} y {1} ".format(px, py))
        print("number_of_pixels {0}".format(number_of_pixels))
        plt.subplot(121)
        plt.imshow(data, origin="lower", vmax=0, vmin=-10)

    y = np.linspace(0, (imgSize - 1), imgSize)
    x = np.linspace(0, (imgSize - 1), imgSize)
    x, y = np.meshgrid(x, y)

    array_radius = np.sqrt((x - px)**2 + (y-py)**2 )
    array_radius[array_radius > number_of_pixels] = -10
    array_radius[array_radius != -10] = 1
    array_radius[array_radius != 1] = 0

    if plot:
        plt.subplot(122)
        plt.imshow(array_radius, origin="lower")
        plt.show()
    
    return data * array_radius

def twoD_Gaussian(data_tuple, amplitude, xo, yo, sigma):
    """
    model of a 2D circular beam

    Parameters
    ----------
    data_tupe   : the x and y mesh grid
    amplitude   : the amp at the peak
    xo          : the peak x location
    yo          : the peak y location
    sigma       : the variance of circular beam

    Returns
    -------
    returns a linear array of the beam
    """

    (x, y) = data_tuple
    xo = float(xo)
    yo = float(yo)    
    a = (1**2)/(2*sigma**2) 
    b =  0
    c =  (1)/(2*sigma**2)
    g = amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()



def getHeadPeaks(headTimeSteps, noradid, mwa, sat, ts):
    """
    finds all the peak pixels in head using gaussian fit

    Parameters
    ----------
    headTimeSteps   : the list of head time Steps
    noradid         : the norad id of the sat
    mwa             : observer object
    sat             : sat object
    ts              : time object

    Returns
    -------
    x_array         : array of x pixels of head
    y_array         : array of y pixels of head
    ra_array        : array of ra of pixels
    dec_array       : array of dec of pixels
    """
    x_array, y_array, ra_array, dec_array = [], [], [], []

    hduBeam = fits.open(beamFile)
    beam = hduBeam[0].data

  
    if debug:
        eprint("finding head pixels")
    
    y = np.linspace(0, 4, 5)
    x = np.linspace(0, 4, 5)
    x, y = np.meshgrid(x, y)  

    for t in tqdm(headTimeSteps):
        hdu = fits.open(prefix + "-t" + str(t).zfill(4) + ".fits")
        UTCTime = datetime.strptime(hdu[0].header['DATE-OBS'][:-2], '%Y-%m-%dT%H:%M:%S')
        data = hdu[0].data
        data = maskData(data, UTCTime, mwa, sat, wcs, ts)

        if np.all(data == 0):
            eprint("file full of zeros")
            eprint("aborting...")
            sys.exit(0)

        row, col = np.where(data == np.max(data))
        beam_cutout = Cutout2D(beam, (col, row), (5 ,5 ))
        cutout = Cutout2D(data, (col, row), (5,5))
        temp = cutout.data / beam_cutout.data
        temp /= np.sum(temp)

        initial_guess = (temp[1,1], 1, 1, 2)
        popt, pconv = opt.curve_fit(twoD_Gaussian, (x, y), temp.ravel(), p0=initial_guess)
        perr = np.sqrt(np.diag(pconv))
        deltaX, deltaY = popt[1], popt[2]
        solX = col + (deltaX - 2)
        solY = row + (deltaY - 2)

        local_x = ufloat(solX[0], perr[1])
        local_y = ufloat(solY[0], perr[2])
        x_array.append(local_x)
        y_array.append(local_y)
        pixcrd = np.array([[0 ,0 ],[solX, solY]], dtype=np.float64)
        world = wcs.wcs_pix2world(pixcrd, 0)
        ra, dec = world[1]
        ra_array.append(ra)
        dec_array.append(dec)

    return x_array, y_array, ra_array, dec_array


def getTailPeaks(TailTimeSteps, noradid, mwa, sat, ts):
    """
    finds all the peak pixels in tail using gaussian fit

    Parameters
    ----------
    TeadTimeSteps   : the list of tail time Steps
    noradid         : the norad id of the sat
    mwa             : observer object
    sat             : sat object
    ts              : time object

    Returns
    -------
    x_array         : array of x pixels of tail
    y_array         : array of y pixels of tail
    ra_array        : array of ra of pixels
    dec_array       : array of dec of pixels
    """

    x_array, y_array, ra_array, dec_array = [], [], [], []
    hduBeam = fits.open(beamFile)
    beam = hduBeam[0].data

    if debug:
        eprint("finding Tail pixels")

    y = np.linspace(0, 4, 5)
    x = np.linspace(0, 4, 5)
    x, y = np.meshgrid(x, y)

    for t in tqdm(TailTimeSteps):

        hdu = fits.open("Neg" + prefix + "-t" + str(t).zfill(4) + ".fits")
        UTCTime = datetime.strptime(hdu[0].header['DATE-OBS'][:-2], '%Y-%m-%dT%H:%M:%S')
        data = hdu[0].data 
        data = maskData(data, UTCTime, mwa, sat, wcs, ts)

        if np.all(data == 0):
            eprint("file full of zeros")
            tailTimeSteps.remove(t)
            continue


        row, col = np.where(data == np.min(data))
        beam_cutout = Cutout2D(beam, (col, row), (5 ,5 ))
        cutout = Cutout2D(data, (col, row), (5,5))
        temp = cutout.data / beam_cutout.data
        temp /= abs(np.sum(temp))
        
        initial_guess = (temp[1,1], 1, 1, 2)
        popt, pconv = opt.curve_fit(twoD_Gaussian, (x, y), temp.ravel(), p0=initial_guess)
        perr = np.sqrt(np.diag(pconv))
        deltaX, deltaY = popt[1], popt[2]
        solX = col + (deltaX - 2)
        solY = row + (deltaY - 2)

        local_x = ufloat(solX[0], perr[1])
        local_y = ufloat(solY[0], perr[2])
        x_array.append(local_x)
        y_array.append(local_y)

        pixcrd = np.array([[0 ,0 ],[solX, solY]], dtype=np.float64)
        world = wcs.wcs_pix2world(pixcrd, 0)
        ra, dec = world[1]
        ra_array.append(ra)
        dec_array.append(dec)

    return x_array, y_array, ra_array, dec_array


def fitLine(head_x, head_y, tail_x, tail_y,order=2):
    
    combined_x_u = head_x + tail_x
    combined_y_u = head_y + tail_y

    combined_x = []
    combined_y = []
    x_errors = []
    y_errors = []

    def poly(x, a, b, c):
        return a*x**2 + b*x + c

    for x, y in zip(combined_x_u, combined_y_u):
        combined_x.append(x.nominal_value)
        combined_y.append(y.nominal_value)
        x_errors.append(x.std_dev)
        y_errors.append(y.std_dev)


    if motion == "E-W":
        # y = f(x)
        ## check if x is going min-max or max-min

        if head_x[0] - head_x[-1] > 0:
            # motion x max to min
            # sort the values max to min
            combined_y = [x for _,x in sorted(zip(combined_x, combined_y),reverse=True)]
            combined_x.sort(reverse=True)
            
        elif head_x[0] - head_x[-1] < 0:
            # motion x min to max
            # sort the values min to max
            combined_y = [x for _,x in sorted(zip(combined_x, combined_y))]
            combined_x.sort()
        
        else:
            eprint("bug")
            eprint("cannot classify motion as min to max for E-W")
            eprint("aborting...")
            sys.exit(0)

        x_array = np.linspace(combined_x[0], combined_x[-1], 100)
        combined_x = np.array(combined_x)
        combined_y = np.array(combined_y)
        z = np.polyfit(combined_x.ravel(), combined_y.ravel(), order)
        f = np.poly1d(z)
        output = [x_array, f(x_array), f]
    
    if motion == "N-S":
        # x = f(y)
        ## check if y is going min-max or max-min

        if head_y[0] - head_y[-1] > 0:
            #motion y max to min
            # sort the values max to min
            combined_x = [x for _,x in sorted(zip(combined_y, combined_x), reverse=True)]
            combined_y.sort(reverse=True)

        elif head_y[0] - head_y[-1] < 0:
            # motion y min to max
            # sort values min to max
            combined_x = [x for _,x in sorted(zip(combined_y, combined_x))]
            combined_y.sort()

        else:
            eprint("bug")
            eprint("cannot classify motion as min to max for N-W")
            eprint("aborting")
            sys.exit(0)

        y_array = np.linspace(combined_y[0], combined_y[-1], 100)
        combined_x = np.array(combined_x)
        combined_y = np.array(combined_y)


        z = np.polyfit(combined_y.ravel(), combined_x.ravel(), order)
        print("z")
        print(z)
        f = np.poly1d(z)
        print("f")
        print(f)
        output = [f(y_array), y_array, f]

    return output



def main(obs, norad):

    # get relevant timeSteps
    global headTimeSteps, tailTimeSteps, midTimeSteps
    headTimeSteps, tailTimeSteps, midTimeSteps = obtainTimeSteps(norad)

    if debug:
        eprint("the selected Headtimesteps " + str(headTimeSteps))
        eprint("the selected Tailtimesteps " + str(tailTimeSteps))
        eprint("the selected Midtimesteps " + str(midTimeSteps))
    
    # get global variables
    global wcs, imgSize, refUTC, imgScale
    hdu = fits.open(prefix + "-t" + str(headTimeSteps[0]).zfill(4) + ".fits")
    wcs = WCS(hdu[0].header, naxis=2)
    imgSize = hdu[0].header["NAXIS1"]
    imgScale = np.abs(hdu[0].header["CDELT2"])
    refUTC = datetime.strptime(hdu[0].header['DATE-OBS'][:-2], '%Y-%m-%dT%H:%M:%S')

    # obtain tle for the date range
    line1, line2, line3 = obtainTLE(norad)


    #get satellite object
    sat, mwa, ts = getSatMWA(line1, line2, line3)    

    if debug:
        eprint("line 1 {0}\nline 2 {1}\nline 3 {2}".format(line1, line2, line3))  

    ## obtain the head data
    head_x, head_y, head_ra, head_dec = getHeadPeaks(headTimeSteps, norad, mwa, sat, ts)
    tail_x, tail_y, tail_ra, tail_dec = getTailPeaks(tailTimeSteps, norad, mwa, sat, ts)


    ## determine if E-W or N-S motion
    deltaX = abs(max(head_x) - min(head_x))
    deltaY = abs(max(head_y) - min(head_y))
    global motion
    if deltaX >= deltaY:
        motion = "E-W"
    elif deltaY > deltaX:
        motion = "N-S"

    if debug:
        eprint("motion classified as " + motion)

    # fit line through head data and extract the mid points
    x_fit, y_fit , function = fitLine(head_x, head_y, tail_x, tail_y)


if __name__ == "__main__":
    parser = ArgumentParser("OrbInfo", description="extracts all the info for orbit determination")
    parser.add_argument("--obs", required=True, type=int, help="the observation id")
    parser.add_argument("--noradid",required=True, type=int, help="The norad id")
    parser.add_argument("--beam", required=True, help="the beam file")
    parser.add_argument("--user", required=True, help="User name for space-track.org")
    parser.add_argument("--passwd", required=True, help="Password for space-track.org")
    parser.add_argument("--debug",default=False,type=bool,help="run script in debug mode")
    parser.add_argument("--filePrefix", default="6Sigma1FloodfillSigmaRFIBinaryMapPeakFlux", type=str, help="the prefix of the file name")
    args = parser.parse_args()

    global debug, query, beamFile, prefix
    debug = args.debug
    query = st.SpaceTrackClient(args.user, args.passwd)
    prefix = args.filePrefix
    beamFile = args.beam

    if debug:
        eprint("running in debug mode")

    main(args.obs, args.noradid)


    
