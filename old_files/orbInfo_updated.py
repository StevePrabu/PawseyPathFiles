#!/usr/bin/env python
from __future__ import division, print_function
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
import matplotlib.pyplot as plt
from astropy.coordinates import EarthLocation
from astropy.time import Time
from astropy import units as u
import csv
import json
from os import path
import sys


def obtainTimeSteps(args):
    with open("filtered_summary.json") as f:
        filtered_summary = json.load(f)

    for line in filtered_summary:
        if filtered_summary[line]["norad"] == str(args.noradid):
            timeSteps = list(map(int,filtered_summary[line]["timeSteps"].split("-")))
    
    headTimeSteps, tailTimeSteps, midTimeSteps = [], [], []
    for t in timeSteps:
        headTimeSteps.append(t)
        tailTimeSteps.append(t+1)

    midTimeSteps = [i for i in headTimeSteps if i in tailTimeSteps]

    return headTimeSteps, tailTimeSteps, midTimeSteps



def obtainTLE(noradid):

    time2 = refUTC + timedelta(hours=24)
    day1, month1, year1 = str(refUTC.day).zfill(2), str(refUTC.month).zfill(2), str(refUTC.year)
    day2, month2, year2 = str(time2.day).zfill(2), str(time2.month).zfill(2), str(time2.year)
    date_range = ops.make_range_string(year1 + "-" + month1 + "-" + day1, year2 + "-" + month2 + "-" + day2)
    debug = True
    if path.exists(str(noradid) + ".txt") == False:

        if debug:
            print("requesting file from server")

        result = query.tle_query(epoch=date_range,norad_cat_id=noradid)

        with open(str(noradid) + ".txt", "w") as outfile:
            json.dump(result.json(), outfile)

        line1 = result.json()[0]["OBJECT_NAME"]
        line2 = result.json()[0]["TLE_LINE1"]
        line3 = result.json()[0]["TLE_LINE2"]

    else:

        if debug:
            print("tle file found on disk. Not downloading.")

        with open(str(noradid) + ".txt") as json_file:
            result = json.load(json_file)

        line1 = result[0]["OBJECT_NAME"]
        line2 = result[0]["TLE_LINE1"]
        line3 = result[0]["TLE_LINE2"]


    return line1 , line2, line3


def getSatMWA(line1, line2, line3):
    ts = load.timescale()
    satellite = EarthSatellite(line2, line3, line1, ts)
    mwa = Topos("26.701276778 S", "116.670846137 E", elevation_m= 377.827)
    return satellite, mwa, ts


def maskData(data, UT, mwa, sat, wcs, ts):
    time = ts.utc(UT.year, UT.month,  UT.day, UT.hour, UT.minute, UT.second)
    sat.at(time)
    difference = sat - mwa
    topocentric = difference.at(time)
    ra, dec, distance = topocentric.radec()
    ra = np.degrees(ra.radians)
    dec = np.degrees(dec.radians)
    px, py = wcs.all_world2pix([[ra, dec]], 1)[0]
    LOS_range = distance.m

    radius = np.degrees(70000/LOS_range) # cone radius (s = r * theta)
    number_of_pixels = radius/imgScale


    y = np.linspace(0, (imgSize - 1), imgSize)
    x = np.linspace(0, (imgSize - 1), imgSize)
    x, y = np.meshgrid(x, y)
    array_radius = np.sqrt((x - px)**2 + (y-py)**2 )
    array_radius[array_radius > number_of_pixels] = -10
    array_radius[array_radius != -10] = 1
    array_radius[array_radius != 1] = 0
    return data * array_radius

def twoD_Gaussian(data_tuple, amplitude, xo, yo, sigma):

    (x, y) = data_tuple
    xo = float(xo)
    yo = float(yo)    
    a = (1**2)/(2*sigma**2) 
    b =  0
    c =  (1)/(2*sigma**2)
    g = amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()


def getTailPeaks(tailTimeSteps, noradid, mwa, sat, ts):
    x_array, y_array, ra_array, dec_array = [], [], [], []

    hduBeam = fits.open(beamFile)
    beam = hduBeam[0].data

    y = np.linspace(0, 4, 5)
    x = np.linspace(0, 4, 5)
    x,y = np.meshgrid(x, y)

    for t in tailTimeSteps:
        hdu = fits.open("Neg" + "6Sigma1Floodfilllr14SigmaRFIBinaryMapPeakFlux-t"+str(t).zfill(4)+ ".fits")
        UTCTime = datetime.strptime(hdu[0].header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
        data = hdu[0].data
        data = maskData(data, UTCTime, mwa, sat, wcs, ts)
        if np.all(data == 0):
             tailTimeSteps.remove(t)
             continue

        row, col = np.where(data == np.min(data))
        beam_cutout = Cutout2D(beam, (col, row), (5, 5))
        cutout = Cutout2D(data, (col, row), (5, 5))
        temp = cutout.data / beam_cutout.data
        temp /= abs(np.sum(temp))

        initial_guess = (temp[1,1], 1, 1, 2)
        popt, pconv = opt.curve_fit(twoD_Gaussian, (x, y), temp.ravel(), p0=initial_guess)
        perr = np.sqrt(np.diag(pconv))
        mu_x, mu_y = popt[1], popt[2]
        print("tail t {} x {} y {} dx {} dy {}".format(t, mu_x, mu_y, perr[1], perr[2]))

        solX = col + (mu_x - 2)
        solY = row + (mu_y - 2)

        x_array.append(solX[0])
        y_array.append(solY[0])
        pixcrd = np.array([[0 ,0 ],[solX, solY]], dtype=np.float64)
        world = wcs.wcs_pix2world(pixcrd,0)
        ra, dec = world[1]
        ra_array.append(ra)
        dec_array.append(dec)

    return x_array, y_array, ra_array, dec_array


def fitLine(head_x, head_y, tail_x, tail_y):
    order = 2
    combined_x = head_x + tail_x
    combined_y = head_y + tail_y

    if motion == "E-W":
        if head_x[0] - head_x[-1] > 0:
            combined_y = [x for _,x in sorted(zip(combined_x, combined_y),reverse=True)]  
            combined_x.sort(reverse=True)

        elif head_x[0] - head_x[-1] < 0:
            combined_y = [x for _,x in sorted(zip(combined_x, combined_y))]
            combined_x.sort()

        else:
            print("BUG")
            print("exiting..")
            sys.exit(0)

        x_array = np.linspace(combined_x[0], combined_x[-1], 100)
        combined_x = np.array(combined_x)
        combined_y = np.array(combined_y)
        z = np.polyfit(combined_x.ravel(), combined_y.ravel(), order)
        f = np.poly1d(z)
        output = [x_array, f(x_array), f]


    if motion == "N-S":
        if head_y[0] - head_y[-1] > 0:
            combined_x = [x for _,x in sorted(zip(combined_y, combined_x), reverse=True)]
            combined_y.sort(reverse=True)
        elif head_y[0] - head_y[-1] < 0:
            combined_x = [x for _,x in sorted(zip(combined_y, combined_x))]
            combined_y.sort()

        else:
            print("BUG")
            print("exiting..")
            sys.exit(0)
        y_array = np.linspace(combined_y[0], combined_y[-1], 100)
        combined_x = np.array(combined_x)
        combined_y = np.array(combined_y)
        z = np.polyfit(combined_y.ravel(), combined_x.ravel(), order)
        f = np.poly1d(z)
        output = [f(y_array), y_array, f]
    return output




def getHeadPeaks(headTimeSteps, noradid, mwa, sat, ts):

    x_array, y_array, ra_array, dec_array = [], [], [], []

    hduBeam = fits.open(beamFile)
    beam = hduBeam[0].data


    y = np.linspace(0, 4, 5)
    x = np.linspace(0, 4, 5)
    x, y = np.meshgrid(x, y)

    for t in headTimeSteps:
        hdu = fits.open("6Sigma1Floodfilllr14SigmaRFIBinaryMapPeakFlux-t"+str(t).zfill(4)+ ".fits")
        UTCTime = datetime.strptime(hdu[0].header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
        data = hdu[0].data
        data = maskData(data, UTCTime, mwa, sat, wcs, ts)
        row, col = np.where(data == np.max(data))
        beam_cutout = Cutout2D(beam, (col, row), (5, 5))
        cutout = Cutout2D(data, (col, row), (5, 5))
        temp = cutout.data/ beam_cutout.data
        temp/= np.sum(temp)

        initial_guess = (temp[1,1], 1, 1, 2)
        popt, pconv = opt.curve_fit(twoD_Gaussian, (x,y), temp.ravel(), p0=initial_guess)
        perr = np.sqrt(np.diag(pconv))
        mu_x, mu_y = popt[1], popt[2]
        print("head t {} x {} y {} dx {} dy {}".format(t, mu_x, mu_y, perr[1], perr[2]))

        solX = col + (mu_x - 2)
        solY = row + (mu_y - 2)

        x_array.append(solX[0])
        y_array.append(solY[0])
        pixcrd = np.array([[0 ,0 ],[solX, solY]], dtype=np.float64)
        world = wcs.wcs_pix2world(pixcrd,0)
        ra, dec = world[1]
        ra_array.append(ra)
        dec_array.append(dec)

    return x_array, y_array, ra_array, dec_array

def s2UTC(s):
    return refUTC + timedelta(seconds=s)


def UTC2s(time):
    return (time - refUTC).total_seconds()


def getMidPoints(midTimeSteps, noradid, mwa, sat, ts, x_fit, y_fit, function):
    x_array, y_array, ra_array, dec_array = [],[], [], []
    time_array = []
    ## update mid points
    midTimeSteps = [i for i in headTimeSteps if i in tailTimeSteps]
    print("updated midTimeSteps {}".format(midTimeSteps))
 
    for t in midTimeSteps:
        hdu_head = fits.open("6Sigma1Floodfilllr14SigmaRFIBinaryMapPeakFlux-t"+str(t).zfill(4)+ ".fits")
        hdu_tail = fits.open("Neg" + "6Sigma1Floodfilllr14SigmaRFIBinaryMapPeakFlux-t"+str(t).zfill(4)+ ".fits")

        head_data = hdu_head[0].data
        head_data /= np.abs(np.sum(head_data))
        tail_data = hdu_tail[0].data
        tail_data /= np.abs(np.sum(tail_data))

        UTCTime = datetime.strptime(hdu_head[0].header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
        time_array.append(UTC2s(UTCTime))

        head_data = maskData(head_data, UTCTime, mwa, sat, wcs, ts)
        tail_data = maskData(tail_data, UTCTime, mwa, sat, wcs, ts)

        streak = head_data + tail_data

        ## iterate through fit and fine boundarys
        fit_val = []
        x_inv, y_inv = [] ,[]
        for x, y in zip(x_fit, y_fit):
            fit_val.append(streak[int(y), int(x)])

        if motion == "E-W":
            min_arg = np.where(fit_val == np.min(fit_val))[0]
            max_arg = np.where(fit_val ==  np.max(fit_val))[0]

            x_min = x_fit[min_arg]
            x_max = x_fit[max_arg]
       
            x_streak = np.linspace(x_min, x_max, 100)
            y_streak = function(x_streak)


            fit_val = []
            for x,y in zip(x_streak, y_streak):
                fit_val.append(streak[int(y), int(x)])

            mid_arg = np.where(np.abs(fit_val) ==  np.min(np.abs(fit_val)))[0]
            if len(mid_arg) == 1:
                x_mid = x_streak[mid_arg]
                y_mid = y_streak[mid_arg]
            elif len(mid_arg) > 1:
                x_mid = x_streak[mid_arg[0]]
                y_mid = y_streak[mid_arg[0]]

        if motion == "N-S":
            min_arg = np.where(fit_val == np.min(fit_val))[0]
            max_arg = np.where(fit_val == np.max(fit_val))[0]

            if float(min_arg.shape[0]) > 1:
                min_arg = min_arg[0]
            elif float(max_arg.shape[0]) > 1:
                max_arg = max_arg[0]

            y_min = y_fit[min_arg]
            y_max = y_fit[max_arg]


            y_streak = np.linspace(y_min, y_max, 100)
            x_streak = function(y_streak)
    
            fit_val = []
            for x,y in zip(x_streak, y_streak):
                fit_val.append(streak[int(y), int(x)])

            mid_arg = np.where(np.abs(fit_val) == np.min(np.abs(fit_val)))[0]
            if len(mid_arg) == 1:
                x_mid = x_streak[mid_arg]
                y_mid = y_streak[mid_arg]
            elif len(mid_arg) > 1:
                x_mid = x_streak[mid_arg[0]]
                y_mid = y_streak[mid_arg[0]]

        x_mid = float(x_mid)
        y_mid = float(y_mid)
        x_array.append(x_mid)
        y_array.append(y_mid)
        pixcrd = np.array([[x_mid,y_mid]],dtype=np.float64)
        world = wcs.wcs_pix2world(pixcrd, 0)[0]
        ra, dec = world
        ra_array.append(ra)
        dec_array.append(dec)
        print("mid t {} x {} y {}".format(t, x_mid, y_mid))

    return x_array, y_array, ra_array, dec_array, time_array





def getUTC(mid_time):
    mid_UTC = []
    for mt in mid_time:
        mid_UTC.append(s2UTC(mt))
    return mid_UTC


def main(args):
    ## get relevant timeSteps
    global headTimeSteps, tailTimeSteps, midTimeSteps
    headTimeSteps, tailTimeSteps, midTimeSteps = obtainTimeSteps(args)

    print("the selected headTimeSteps " + str(headTimeSteps))
    print("the selected tailTimeSteps " + str(tailTimeSteps))
    print("the selected midTimeSteps " + str(midTimeSteps))

    ## get global variables
    global wcs, imgSize, refUTC, imgScale
    hdu = fits.open("6Sigma1Floodfilllr14SigmaRFIBinaryMapSeed-t"+ str(headTimeSteps[0]).zfill(4) + ".fits")
    wcs = WCS(hdu[0].header, naxis=2)
    imgSize = hdu[0].header["NAXIS1"]
    imgScale = np.abs(hdu[0].header["CDELT2"])
    refUTC = datetime.strptime(hdu[0].header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')

    ## obtain tle for the date range
    line1, line2, line3 = obtainTLE(args.noradid)

    # get satellite objeect
    sat, mwa, ts = getSatMWA(line1, line2, line3)

    ##  obtain head data
    head_x, head_y, head_ra, head_dec = getHeadPeaks(headTimeSteps, args.noradid, mwa, sat, ts)    

    ## obtain tail data
    tail_x, tail_y, tail_ra, tail_dec = getTailPeaks(tailTimeSteps, args.noradid, mwa, sat, ts)

    ## determine if motion E-W or N-S
    deltaX = abs(np.max(head_x) - np.min(head_x))
    deltaY = abs(np.max(head_y) - np.min(head_y))
    global motion
    if deltaX >= deltaY:
        motion = "E-W"
    elif deltaY > deltaX:
        motion = "N-S"

    print("motion classified as {}".format(motion))
    ## fit line through head data and extrac the mid points
    x_fit, y_fit, function = fitLine(head_x, head_y, tail_x, tail_y)
    
    ### get mid points
    mid_x, mid_y, mid_ra, mid_dec, mid_time = getMidPoints(midTimeSteps, args.noradid, mwa, sat, ts, x_fit, y_fit, function) 

    mid_x = np.array(mid_x)
    mid_y = np.array(mid_y)
    mid_ra = np.array(mid_ra)
    mid_dec = np.array(mid_dec)
    mid_time = np.array(mid_time)

    ## save values
    np.save(str(args.noradid) + "_x_array.npy", mid_x)
    np.save(str(args.noradid) + "_y_array.npy", mid_y)
    np.save(str(args.noradid) + "_ra_array.npy", mid_ra)
    np.save(str(args.noradid) + "_dec_array.npy", mid_dec)


    ## save extracted info to csv file
    mid_UTC = getUTC(mid_time)
    with open(str(args.noradid) + "_extracted_data_from_" + str(args.obs) + ".csv", "w") as vsc:
        thewriter = csv.writer(vsc)
        thewriter.writerow(["x", "y", "ra", "dec" , "UTC", "t" ])
        for x, y, ra, dec, utc, t in zip(mid_x, mid_y, mid_ra, mid_dec, mid_UTC, mid_time):
            output = [x, y, ra, dec, utc, t]

            thewriter.writerow(output)


if __name__ == "__main__":
    parser = ArgumentParser("OrbInfoUpdated", description="extracts all the info for orbit determination")
    parser.add_argument("--obs",required=True,type=int,help="the observatio id")
    parser.add_argument("--noradid",required=True,type=int,help="the noradid of the satellite")
    parser.add_argument("--beam",required=True,help="the beam file")
    parser.add_argument("--user",required=True, help="the user name for spacetrack.org")
    parser.add_argument("--passwd",required=True,help="the passord for spacetrack.org")
    args = parser.parse_args()
    global beamFile, query
    query = st.SpaceTrackClient(args.user, args.passwd)
    beamFile = args.beam

    main(args)
