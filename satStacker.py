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
import matplotlib.pyplot as plt
from os import path
import json
from scipy.fft import fft
from scipy import signal as s
from multiprocessing import Pool
import multiprocessing

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


def getXY(sat, mwa, utc_ts):

    sat.at(utc_ts)
    difference = sat - mwa
    topocentric = difference.at(utc_ts)
    ra, dec, distance = topocentric.radec()
    ra = np.degrees(ra.radians)
    dec = np.degrees(dec.radians)
    pix_coords = wcs.all_world2pix(np.array([[ra], [dec]]).T, 1)

    return pix_coords.T

def getsnr(diff, x, y, window):
    cutout = Cutout2D(diff, (x,y), (window, window))
    
    data = cutout.data

    ## do 2 rounds of source masking to calculate accurate local rms
    rms1 = np.std(data)
    data[np.abs(data) > 3*rms1] = 0
    rms2 = np.std(data)
    data[np.abs(data) > 3*rms2] = 0

    rms = np.std(data)
  
    return diff[int(y),int(x)]/rms


def lowpass(signal):
    T = 2
    N = len(signal)
    fc = 0.1 ## cutoff freq for 10s long signal
    w = fc/(0.5/2)
    b,a = s.butter(5, w, "low")
    return s.filtfilt(b,a,signal)

def integrate(signal):
    signal = np.array(signal)
    norm_factor = np.count_nonzero(signal)
    if norm_factor == 0:
        norm_factor = 1
    filtered_data = lowpass(signal)

    return np.sum(filtered_data)/norm_factor


#def worker(phase):
#    signal = []
#    for t in range(args.timeSteps):
#        hdu1 = fits.open(str(args.obs) + "-2m-"+str(t)+"-"+str(f).zfill(4)+"-dirty.fits")
#        hdu2 = fits.open(str(args.obs) + "-2m-"+str(t+1)+"-"+str(f).zfill(4)+"-dirty.fits")
#        diff = hdu2[0].data[0,0,:,:] - hdu1[0].data[0,0,:,:]
#        
#        if np.all(diff == 0):
#            signal.append(0)
#            continue
#        
#        utc = datetime.strptime(hdu2[0].header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f') + timedelta(seconds=phase)
#        hdu1.close()
#        hdu2.close()
#        utc_ts = ts.utc(utc.year, utc.month, utc.day, utc.hour, utc.minute, utc.second)
#        x, y = getXY(sat, mwa, utc_ts)
#        if np.sqrt((x-imgSize/2)**2 + (y-imgSize/2)**2)*imgScale > 18:
#            continue
#        #snr = getsnr(diff, x, y, args.window)
#        snr = diff[int(y), int(x)]
#        signal.append(snr)
#    output = integrate(signal)
#    return output

def worker(phase):
    
    signal = []
    counter = 0
    for utc in utc_array:
        local_utc = utc + timedelta(seconds=phase)
        utc_ts = ts.utc(local_utc.year, local_utc.month, local_utc.day, local_utc.hour, local_utc.minute, local_utc.second)
        x, y = getXY(sat, mwa, utc_ts)
        if np.sqrt((x-imgSize/2)**2 + (y-imgSize/2)**2)*imgScale > 18:
            signal.append(0)
            continue
        snr = cube[counter,int(y),int(x)]
        signal.append(snr)
        counter += 1
    
    output = integrate(signal)
    return output


def getCube(f,args):
    if debug:
        print("creting cube")
    output_cube = np.zeros((args.timeSteps,imgSize,imgSize))
    output_utc_array = []
    for t in range(args.timeSteps):
        if debug:
            print("reading t {}".format(t))
        hdu1 = fits.open(str(args.obs) + "-"+str(args.midName)+"-"+str(t)+"-"+str(f).zfill(4)+"-dirty.fits")
        hdu2 = fits.open(str(args.obs) + "-"+str(args.midName)+"-"+str(t+1)+"-"+str(f).zfill(4)+"-dirty.fits")
        diff = hdu2[0].data[0,0,:,:] - hdu1[0].data[0,0,:,:]
        utc = datetime.strptime(hdu2[0].header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
        hdu1.close()
        hdu2.close()
        output_utc_array.append(utc)
        output_cube[t,:,:] = diff

    if debug:
        print("done")
    return output_utc_array, output_cube

def main(args):
    ## get the global parameters
    global wcs, imgSize, refUTC, imgScale, sat, mwa, ts
    hdu = fits.open(str(args.obs) + "-"+str(args.midName)+"-1-0000-dirty.fits" )
    wcs = WCS(hdu[0].header, naxis=2)
    imgSize = hdu[0].header["NAXIS1"]
    imgScale = np.abs(hdu[0].header["CDELT2"])
    refUTC = datetime.strptime(hdu[0].header["DATE-OBS"], '%Y-%m-%dT%H:%M:%S.%f')

    ## get the relevant tle for the satellite
    line1, line2, line3 = obtainTLE(args.noradid)

    if debug:
        print("line1 {0}\nline2 {1}\nline3 {2}".format(line1, line2, line3))


    ## get the satellite and observer object
    sat, mwa, ts = getSatMWA(line1,  line2, line3)

    phase_array = np.linspace(-10,10,20)

    waterfall = np.zeros((20,args.channels))
    global f, utc_array, cube
    for f in tqdm(range(args.channels)):
        utc_array, cube = getCube(f,args)
        with multiprocessing.Pool(processes=20) as pool:
            results = pool.starmap(worker, np.array([phase_array]).T)
        waterfall[:,f] = results

    np.save(str(args.noradid) + "-" +str(args.obs) + "waterfall.npy", waterfall)
    plt.imshow(waterfall, origin="lower", aspect="auto")
    plt.show()

     

    

                 
if __name__ == "__main__":
    global args
    parser = ArgumentParser("satStack", description="stacks and intergrates over the snr of sat")
    parser.add_argument("--obs", required=True, type=int, help="the observation id")
    parser.add_argument("--noradid", required=True, type=int, help="the noradid of the satellite")
    parser.add_argument("--window",type=int, default=100, help="the box to calculate local noise")
    parser.add_argument("--debug",type=bool,default=False,help="run script in debug mode")
    parser.add_argument("--timeSteps",default=54,type=int, help="the number of timeSteps")
    parser.add_argument("--channels",default=768,type=int,help="the number of fine channels")
    parser.add_argument("--user",required=True, help="the user name for spacetrack.org")
    parser.add_argument("--passwd",required=True, help="the password for spacetrack.org")
    parser.add_argument("--midName",default="2m", help="the mid name of the fits file")
    args = parser.parse_args()
    
    global debug, query
    debug = args.debug
    query = st.SpaceTrackClient(args.user, args.passwd)

    if debug:
        print("running in debug mode")

    main(args)
