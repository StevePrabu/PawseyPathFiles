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
#import matplotlib
#matplotlib.use('Agg')
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

obs_loc = EarthLocation(lon=116.67083333*u.deg, lat=-26.70331941*u.deg, height=377.827*u.m)



def mc_sampler(mu, sigma):
    
    s = np.random.normal(mu, sigma, 200)
    s =  [x for x in s if x > (mu-2*sigma) and x <= (mu+2*sigma)]
    output = float(np.random.choice(s,1)[0])
    if actual:
        output = mu
    weight = 1/sigma
    
    return output, weight

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

    ts = load.timescale()
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
    radius = np.degrees(70000/LOS_range) # cone radius (s = r * theta)
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
    x_weights, y_weights = [], []

    hduBeam = fits.open(beamFile)
    beam = hduBeam[0].data

  
    #if debug:
    #    eprint("finding head pixels")
    
    y = np.linspace(0, 4, 5)
    x = np.linspace(0, 4, 5)
    x, y = np.meshgrid(x, y)  

    for t in headTimeSteps:
        hdu = fits.open(prefix + "-t" + str(t).zfill(4) + ".fits")
        UTCTime = datetime.strptime(hdu[0].header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
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

        mu_x, mu_y = popt[1], popt[2]
        sigma_x, sigma_y = perr[1], perr[2]
        deltaX, x_weight = mc_sampler(mu_x, sigma_x) 
        deltaY, y_weight = mc_sampler(mu_y, sigma_y)
        x_weights.append(x_weight)
        y_weights.append(y_weight)

        solX = col + (deltaX - 2)
        solY = row + (deltaY - 2)
        x_array.append(solX[0])
        y_array.append(solY[0])
        pixcrd = np.array([[0 ,0 ],[solX, solY]], dtype=np.float64)
        world = wcs.wcs_pix2world(pixcrd, 0)
        ra, dec = world[1]
        ra_array.append(ra)
        dec_array.append(dec)

    return x_array, y_array, ra_array, dec_array, x_weights, y_weights

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
    x_weights, y_weights = [], []
    hduBeam = fits.open(beamFile)
    beam = hduBeam[0].data

    #if debug:
    #    eprint("finding Tail pixels")

    y = np.linspace(0, 4, 5)
    x = np.linspace(0, 4, 5)
    x, y = np.meshgrid(x, y)

    for t in TailTimeSteps:

        hdu = fits.open("Neg" + prefix + "-t" + str(t).zfill(4) + ".fits")
        UTCTime = datetime.strptime(hdu[0].header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
        data = hdu[0].data 
        #print(data.shape)
        data = maskData(data, UTCTime, mwa, sat, wcs, ts)

        if np.all(data == 0):
            #eprint("file full of zeros")
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
        mu_x, mu_y = popt[1], popt[2]
        sigma_x, sigma_y = perr[1], perr[2]
        deltaX,x_weight = mc_sampler(mu_x, sigma_x)
        deltaY, y_weight= mc_sampler(mu_y, sigma_y)

        x_weights.append(x_weight)
        y_weights.append(y_weight)

        solX = col + (deltaX - 2)
        solY = row + (deltaY - 2)
        x_array.append(solX[0])
        y_array.append(solY[0])
        pixcrd = np.array([[0 ,0 ],[solX, solY]], dtype=np.float64)
        world = wcs.wcs_pix2world(pixcrd, 0)
        ra, dec = world[1]
        ra_array.append(ra)
        dec_array.append(dec)

    return x_array, y_array, ra_array, dec_array, x_weights, y_weights

def fitLine(head_x, head_y, tail_x, tail_y,order=2):

    combined_x = head_x + tail_x
    combined_y = head_y + tail_y
    w_x = head_x_w + tail_x_w
    w_y = head_y_w + tail_y_w
    #print(w_x)
    #print(w_y)

    if motion == "E-W":
        # y = f(x)
        ## check if x is going min-max or max-min

        if head_x[0] - head_x[-1] > 0:
            # motion x max to min
            # sort the values max to min
            combined_y = [x for _,x in sorted(zip(combined_x, combined_y),reverse=True)]
            w_y = [x for _,x in sorted(zip(combined_x, w_y),reverse=True)]
            combined_x.sort(reverse=True)
            
        elif head_x[0] - head_x[-1] < 0:
            # motion x min to max
            # sort the values min to max
            combined_y = [x for _,x in sorted(zip(combined_x, combined_y))]
            w_y = [x for _,x in sorted(zip(combined_x, w_y))]
            combined_x.sort()
        
        else:
            eprint("bug")
            eprint("cannot classify motion as min to max for E-W")
            eprint("aborting...")
            sys.exit(0)

        x_array = np.linspace(combined_x[0], combined_x[-1], 100)
        combined_x = np.array(combined_x)
        combined_y = np.array(combined_y)
        w_y = np.array(w_y)
        z = np.polyfit(combined_x.ravel(), combined_y.ravel(), order, w=w_y)
        f = np.poly1d(z)
        output = [x_array, f(x_array), f]
    
    if motion == "N-S":
        # x = f(y)
        ## check if y is going min-max or max-min

        if head_y[0] - head_y[-1] > 0:
            #motion y max to min
            # sort the values max to min
            combined_x = [x for _,x in sorted(zip(combined_y, combined_x), reverse=True)]
            w_x = [x for _,x in sorted(zip(combined_y, w_x), reverse=True)]
            combined_y.sort(reverse=True)

        elif head_y[0] - head_y[-1] < 0:
            # motion y min to max
            # sort values min to max
            combined_x = [x for _,x in sorted(zip(combined_y, combined_x))]
            w_x = [x for _,x in sorted(zip(combined_y, w_x))]
            combined_y.sort()

        else:
            eprint("bug")
            eprint("cannot classify motion as min to max for N-W")
            eprint("aborting")
            sys.exit(0)


        y_array = np.linspace(combined_y[0], combined_y[-1], 100)
        combined_x = np.array(combined_x)
        combined_y = np.array(combined_y)
        w_x = np.array(w_x)
        z = np.polyfit(combined_y.ravel(), combined_x.ravel(), order, w=w_x)
        f = np.poly1d(z)
        output = [f(y_array), y_array, f]

    return output

def UTC2s(time):
    """
    converts utc to s with respect to refUTC
    
    Paramters
    ---------
    time    : time to convert into ref seconds

    Returns
    -------
    time converted to ref sec
    """
    
    return (time - refUTC).total_seconds()


def getMidPoints(MidTimeSteps, noradid, mwa, sat, ts, x_fit, y_fit, function):

    x_array, y_array, ra_array, dec_array, time_array = [] ,[], [], [], []

    if debug:
        #eprint("recalculating new mid timeSteps")
        midTimeSteps = [i for i in headTimeSteps if i in tailTimeSteps]
        #eprint("updated timeSteps are ")
        #eprint("head " + str(headTimeSteps))
        #eprint("tail " + str(tailTimeSteps))
        #eprint("mid " + str(midTimeSteps))

    for t in midTimeSteps:
        
        hdu_head = fits.open(prefix+"-t" + str(t).zfill(4) + ".fits")
        hdu_tail = fits.open("Neg" + prefix + "-t" + str(t).zfill(4) + ".fits")

        head_data = hdu_head[0].data
        head_data /= np.abs(np.sum(head_data))
        tail_data = hdu_tail[0].data
        tail_data /= np.abs(np.sum(tail_data))

        UTCTime = datetime.strptime(hdu_head[0].header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')

        time_array.append(UTC2s(UTCTime))

        head_data = maskData(head_data, UTCTime, mwa, sat, wcs, ts)
        tail_data = maskData(tail_data, UTCTime, mwa, sat, wcs, ts)

        streak = head_data + tail_data

        ## iterate throug fit and find boundarys
        fit_val = []
        #x_prev = 
        x_inv, y_inv = [], []
        for x, y in zip(x_fit, y_fit):
            fit_val.append(streak[int(y),int(x)])

        if motion == "E-W":
            min_arg = np.where(fit_val == np.min(fit_val))[0]
            max_arg = np.where(fit_val == np.max(fit_val))[0]
            x_min = x_fit[min_arg]
            x_max = x_fit[max_arg]

            x_streak = np.linspace(x_min, x_max, 100)
            y_streak = function(x_streak)


            fit_val = []
            #x_inv, y_inv = [], []
            for x, y in zip(x_streak, y_streak):
                fit_val.append(streak[int(y), int(x)])
  
            
            mid_arg = np.where(np.abs(fit_val) == np.min(np.abs(fit_val)))[0]
            if len(mid_arg) == 1:
                x_mid = x_streak[mid_arg]
                y_mid = y_streak[mid_arg]
            elif len(mid_arg) > 1:
                x_mid = x_streak[mid_arg[0]]
                y_mid = y_streak[mid_arg[0]]


        if motion == "N-S":
            min_arg = np.where(fit_val == np.min(fit_val))[0]
            max_arg = np.where(fit_val == np.max(fit_val))[0]

            ### check dimension
            if float(min_arg.shape[0]) > 1:
                min_arg = min_arg[0]
            elif float(max_arg.shape[0]) >1:
                max_arg = max_arg[0]

            y_min = y_fit[min_arg]
            y_max = y_fit[max_arg]

            y_streak = np.linspace(y_min, y_max, 100)
            x_streak = function(y_streak)

            fit_val = []
            print(min_arg.shape)
            print("min_arg {0} max_arg {1}".format(min_arg, max_arg))
            for x,y in zip(x_streak, y_streak):
                fit_val.append(streak[int(y), int(x)])

            mid_arg = np.where(np.abs(fit_val) == np.min(np.abs(fit_val)))[0]
            if len(mid_arg) == 1:
                x_mid = x_streak[mid_arg]
                y_mid = y_streak[mid_arg]
            elif len(mid_arg) > 1:
                x_mid = x_streak[mid_arg[0]]
                y_mid = y_streak[mid_arg[0]]

        #print(int(x_mid))
        x_mid = float(x_mid)
        y_mid = float(y_mid)
        x_array.append(x_mid)
        y_array.append(y_mid)
        pixcrd = np.array([[x_mid,y_mid]],dtype=np.float64)
        #pixcrd = np.array([[0 ,0 ],[x_mid, y_mid]], dtype=np.float64)
        world = wcs.wcs_pix2world(pixcrd, 0)[0]
        ra, dec = world
        ra_array.append(ra)
        dec_array.append(dec)

    return x_array, y_array, ra_array, dec_array, time_array

def timeStamp(head_x, head_y, tail_x, tail_y, mid_x, mid_y, mid_time):
    order = 4

    x_z = np.polyfit(mid_x, mid_time, order)
    y_z = np.polyfit(mid_y, mid_time, order)
    x_function = np.poly1d(x_z)
    y_function = np.poly1d(y_z)

    head_time = (x_function(head_x) + y_function(head_y))/2
    tail_time = (x_function(tail_x) + y_function(tail_y))/2

    return head_time, tail_time


def  getUTC(mid_time):
    """
    convert ref time in sec to UTC

    Paramters
    ---------
    head_time   : head ref time array
    tail_time   : tail ref time array
    mid_time    : mid ref time array

    Returns
    -------
    head_UTC    : head utc array
    tail_UTC    : tail utc array
    mid_UTC     : mid utc array
    """

    head_UTC, tail_UTC, mid_UTC  = [], [], []

    #for ht in head_time:
    #    head_UTC.append(s2UTC(ht))

    #for tt in tail_time:
    #    tail_UTC.append(s2UTC(tt))

    for mt in mid_time:
        mid_UTC.append(s2UTC(mt))

    return mid_UTC


def getLST(head_UTC, tail_UTC, mid_UTC):
    """
    obtains lst for the head tail and mid points
    
    Parameters
    ----------
    head_UTC    : array of head UTC
    tail_UTC    : array of tail UTC
    mid_UTC     : array of mid UTC

    Returns
    -------
    head_LST    : array of head lst
    tail_LST    : array of tail lst
    mid_LST     : array of mid lst
    """

    head_LST, tail_LST, mid_LST = [], [], []

    #for ht in head_UTC:
    #    head_LST.append(UTC2LST(ht))

    #for tt in tail_UTC:
    #    tail_LST.append(UTC2LST(tt))

    for mt in mid_UTC:
        mid_LST.append(UTC2LST(mt))
    
    return head_LST, tail_LST, mid_LST


def UTC2LST(UTCTime):
    """
    converts utc to lst angle

    Parameters 
    ----------
    UTCTime     : time to covert to lst

    Returns
    -------
    LST_angle   : lst angle
    """

    obs_time =  Time(UTCTime,location=obs_loc)
    LST = obs_time.sidereal_time("mean")
    h , residue1 = str(LST).split("h")
    m , residue2 = str(residue1).split("m")
    s, residue3 = str(residue2).split("s")
    LST_angle =  15.0*float(h) + 0.25*float(m) + 0.004166666666666667*float(s)
    return LST_angle

def s2UTC(s):
    """
    converts ref s to utc

    Parameters
    ----------
    s   : input ref seconds

    Returns
    -------
    converted utc
    """
    
    return refUTC + timedelta(seconds=s)

def main(obs, norad, niter):

    # get relevant timeSteps
    global headTimeSteps, tailTimeSteps, midTimeSteps, actual
    headTimeSteps, tailTimeSteps, midTimeSteps = obtainTimeSteps(norad)
    actual = False

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
    refUTC = datetime.strptime(hdu[0].header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')

    # obtain tle for the date range
    line1, line2, line3 = obtainTLE(norad)

    #get satellite object
    sat, mwa, ts = getSatMWA(line1, line2, line3)    

    if debug:
        eprint("line 1 {0}\nline 2 {1}\nline 3 {2}".format(line1, line2, line3))  


    ## define all trace elements
    global head_x_w, head_y_w, tail_x_w, tail_y_w
    mid_x_array, mid_y_array, mid_ra_array, mid_dec_array, mid_time_array = [],[],[],[],[]
    for i in tqdm(range(niter)):
        if i == (niter -1):
            print("last iteration")
            actual = True

        ## obtain head data
        head_x, head_y, head_ra, head_dec, head_x_w, head_y_w = getHeadPeaks(headTimeSteps, norad, mwa, sat, ts)
        #print(head_x)
        tail_x, tail_y, tail_ra, tail_dec, tail_x_w, tail_y_w = getTailPeaks(tailTimeSteps, norad, mwa, sat, ts)

        ## determine if E-W or N-S motion
        deltaX = abs(np.max(head_x) - np.min(head_x))
        deltaY = abs(np.max(head_y) - np.min(head_y))
        global motion
        if deltaX >= deltaY:
            motion = "E-W"
        elif deltaY > deltaX:
            motion = "N-S"

        #if debug:
        #    eprint("motion classified as " + motion)

        # fit line through head data and extract the mid points
        x_fit, y_fit , function = fitLine(head_x, head_y, tail_x, tail_y)

        #mid_x, mid_y, mid_ra, mid_dec, mid_time = getMidPoints
        mid_x, mid_y, mid_ra, mid_dec, mid_time = getMidPoints(midTimeSteps, norad, mwa, sat, ts, x_fit, y_fit, function)

        mid_x_array.append(mid_x)
        mid_y_array.append(mid_y)
        mid_ra_array.append(mid_ra)
        mid_dec_array.append(mid_dec)
        mid_time_array.append(mid_time)

        # get timeStamp for the head and tail pixels
        #head_time, tail_time = timeStamp(head_x, head_y, tail_x, tail_y, mid_x, mid_y, mid_time)

        # get utc and lst for data points
        #print("getting utc")

    
    mid_x_array = np.array(mid_x_array)
    mid_y_array = np.array(mid_y_array)
    mid_ra_array = np.array(mid_ra_array)
    mid_dec_array = np.array(mid_dec_array)
    mid_time_array = np.array(mid_time_array)

    ## save trace arrays as numpy array
    np.save(str(norad)+"_x_array.npy", mid_x_array)
    np.save(str(norad)+"_y_array.npy", mid_y_array)
    np.save(str(norad)+"_ra_array.npy", mid_ra_array)
    np.save(str(norad)+"_dec_array.npy", mid_dec_array)
    np.save(str(norad)+"_time_array.npy", mid_time_array)
    ####

    mid_x_val = mid_x_array[-1,:]
    mid_x_max = np.max(mid_x_array, axis=0)
    mid_x_min = np.min(mid_x_array, axis=0)
    
    mid_y_val = mid_y_array[-1,:]
    mid_y_max = np.max(mid_y_array, axis=0)
    mid_y_min = np.min(mid_y_array, axis=0)

    mid_ra_val = mid_ra_array[-1,:]
    mid_ra_max = np.max(mid_ra_array, axis=0)
    mid_ra_min = np.min(mid_ra_array, axis=0)
    
    mid_dec_val = mid_dec_array[-1,:]
    mid_dec_max = np.max(mid_dec_array, axis=0)
    mid_dec_min = np.min(mid_dec_array, axis=0)



    mid_time_val = np.mean(mid_time_array, axis=0)
  
    # get utc
    mid_UTC = getUTC( mid_time_val)
    #head_LST, tail_LST, mid_LST = getLST(head_UTC, tail_UTC, mid_UTC)

    ## save scatter points to file
    with open(str(norad) + "_extracted_mc_data_from_" + str(obs) + ".csv", "w") as vsc:

        if debug:
            eprint("writing the obtain points to file...")        

        thewriter = csv.writer(vsc)
        thewriter.writerow(["x", "y", "ra", "dec", "UTC", "t", "x_max", "y_max", "x_min", "y_min", "ra_max", "ra_min", "dec_max", "dec_min"])

        for x, y, ra, dec, utc, t, x_max, y_max, x_min, y_min, ra_max, ra_min, dec_max, dec_min in zip(mid_x_val, \
        mid_y_val, mid_ra_val, mid_dec_val, mid_UTC, mid_time_val, mid_x_max, mid_y_max, mid_x_min, mid_y_min,\
        mid_ra_max, mid_ra_min, mid_dec_max, mid_dec_min):
            output = [x, y, ra, dec, utc, t, x_max, y_max, x_min, y_min, ra_max, ra_min, dec_max, dec_min]
            thewriter.writerow(output)


    ## save fit points as img
    plt.subplot(111, projection=wcs)
    plt.scatter(head_x, head_y, label="head")
    plt.scatter(tail_x, tail_y, label="tail")
    plt.scatter(mid_x, mid_y, label="mid")
    plt.xlabel("RA (Degrees)")
    plt.ylabel("DEC (Degrees)")
    plt.ylim(0, imgSize)
    plt.xlim(0, imgSize)
    plt.grid()
    plt.legend()
    plt.savefig("scatter_norad"+str(norad)+"_obs" + str(obs)+ ".png", dpi=300)

    #### "That's all Folks!"" ###        
    if debug:
        eprint(" That's all Folks! ")
    

if __name__ == "__main__":
    parser = ArgumentParser("OrbInfo", description="extracts all the info for orbit determination")
    parser.add_argument("--obs", required=True, type=int, help="the observation id")
    parser.add_argument("--noradid",required=True, type=int, help="The norad id")
    parser.add_argument("--beam", required=True, help="the beam file")
    parser.add_argument("--user", required=True, help="User name for space-track.org")
    parser.add_argument("--passwd", required=True, help="Password for space-track.org")
    parser.add_argument("--debug",default=False,type=bool,help="run script in debug mode")
    parser.add_argument("--filePrefix", default="6Sigma1FloodfillSigmaRFIBinaryMapPeakFlux", type=str, help="the prefix of the file name")
    parser.add_argument("--niter",default=100,type=int,help="number of monte carlo iterations")
    args = parser.parse_args()

    global debug, query, beamFile, prefix
    debug = args.debug
    query = st.SpaceTrackClient(args.user, args.passwd)
    prefix = args.filePrefix
    beamFile = args.beam

    if debug:
        eprint("running in debug mode")

    main(args.obs, args.noradid, args.niter)


    
