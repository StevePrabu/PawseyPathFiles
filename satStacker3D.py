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
from multiprocessing import Pool
import multiprocessing

global args

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

        if args.debug:
            print("requesting file from server")

        result = query.tle_query(epoch=date_range,norad_cat_id=noradid)

        with open(str(noradid) + ".txt", "w") as outfile:
            json.dump(result.json(), outfile)

        line1 = result.json()[0]["OBJECT_NAME"]
        line2 = result.json()[0]["TLE_LINE1"]
        line3 = result.json()[0]["TLE_LINE2"]

    else:

        if args.debug:
            print("tle file found on disk. Not downloading.")

        with open(str(noradid) + ".txt") as json_file:
            result = json.load(json_file)

        line1 = result[0]["OBJECT_NAME"]
        line2 = result[0]["TLE_LINE1"]
        line3 = result[0]["TLE_LINE2"]


    return line1 , line2, line3



def constructTLE(i, ra, e, aop, ma, mm, sid):
    
    if ra < 0:
        #print("updating ra {}".format(ra))
        ra += 360
    
    if aop < 0:
        #print("updating aop {}".format(aop))
        aop += 360

    if ma < 0:
        #print("updating ma {}".format(ma))
        ma += 360

    #print(" i {} ra {} e {} aop {} ma {} mm {} ".format(i, ra, e, aop, ma, mm))
    line_template_part1 = "2 " + str(sid).ljust(5) + " "
    line_template_part2 = str("%08.4f" %i) + " " + str("%08.4f" %ra) +\
     " " +  str(str("%09.7f" %e).split(".")[1]) + " " + str("%08.4f" %aop) \
     + " " + str("%08.4f" %ma) + " " + str("%011.8f" %mm)  + "14849"
    line_template_part3 = sumCheck(line_template_part1 + line_template_part2)

    tle = line_template_part1 + line_template_part2 + line_template_part3

    return tle



def sumCheck(digit):

    digit = digit.replace("-", "1")    
    a = sum(int(x) for x in digit if x.isdigit())

    return str(a%10)
    


def updateTLE(lline2,raOffset=0):
    i = float(lline2[8:16])
    ra = float(lline2[17:25]) + raOffset
    e = float(lline2[26:33])/10000000
    aop = float(lline2[34:42])
    ma = float(lline2[43:51])
    mm = float(lline2[52:63])
    return constructTLE(i, ra, e, aop, ma, mm, args.noradid)
    
def getCube(f):

    if args.debug:
        print("reading files to create cube...")
    output_cube = np.zeros((args.timeSteps,imgSize,imgSize))
    
    output_utc = []

    for t in range(args.timeSteps):
        hdu1 = fits.open(str(args.obs) + "-" + str(args.midName) + "-" + str(t) + "-" + str(f).zfill(4) + "-dirty.fits")
        hdu2 = fits.open(str(args.obs) + "-" + str(args.midName) + "-" + str(t+1) + "-" + str(f).zfill(4) + "-dirty.fits")

        data1 = hdu1[0].data[0,0,:,:]
        data2 = hdu2[0].data[0,0,:,:]

        ## check if any of the files are full of zeros
        if np.any(data1==0) or np.any(data2==0):
            diff = np.zeros((1400,1400))
        
        else:
            diff = data2 - data1
       
        
        utc = datetime.strptime(hdu2[0].header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')  
        output_utc.append(utc)
        

        hdu1.close()
        hdu2.close()

        output_cube[t,:,:] = diff

    if args.debug:
        print("done")

    return output_utc, output_cube

def getXY(sat, local_time):
    sat.at(local_time)
    difference = sat - mwa
    topocentric = difference.at(local_time)
    ra, dec, distance = topocentric.radec()
    ra = np.degrees(ra.radians)
    dec = np.degrees(dec.radians)
    pix_coords = wcs.all_world2pix(np.array([[ra], [dec]]).T, 1)

    return pix_coords.T


def worker(intrack_phase, offtrack_phase):

    signal = []
        
    local_line3 = updateTLE(line3,raOffset=offtrack_phase)
    sat = EarthSatellite(line2, local_line3, line1, ts)
    time_counter = 0
    s2_snapshot_time_array = np.linspace(0,2,10)
    x_array, y_array = [], []

    for utc in utc_array:
        
        for deltaT in s2_snapshot_time_array:
            local_utc = utc + timedelta(seconds=(intrack_phase+deltaT))

            utc_ts = ts.utc(local_utc.year, local_utc.month, local_utc.day, local_utc.hour, local_utc.minute, local_utc.second + local_utc.microsecond/1000.0)
            x, y = getXY(sat, utc_ts)

            if np.sqrt((x-imgSize/2)**2 + (y-imgSize/2)**2)*imgScale > 18:
                continue

            if np.isnan(x) or np.isnan(y):                
                continue
           
            signal.append(cube[time_counter, int(y), int(x)])

        time_counter += 1
    if np.all(signal == 0):
        return 0
    elif signal == []:
        return 0
    else:
        return np.mean(signal)


def main():
    ## get the global parameters
    global wcs, imgSize, refUTC, imgScale, mwa, ts, line1, line2, line3
    hdu = fits.open(str(args.obs) + "-" + str(args.midName) + "-1-0000-dirty.fits")
    wcs = WCS(hdu[0].header, naxis=2)
    imgSize = hdu[0].header["NAXIS1"]
    imgScale = np.abs(hdu[0].header["CDELT2"])
    refUTC = datetime.strptime(hdu[0].header["DATE-OBS"], '%Y-%m-%dT%H:%M:%S.%f')

    ## get the relevant tle for the satellite
    line1, line2, line3 = obtainTLE(args.noradid)
    
    if args.debug:
        print("line1 {0}\nline2 {1}\nline3 {2}".format(line1, line2, line3))


    mwa = Topos("26.701276778 S", "116.670846137 E", elevation_m= 377.827)
    ts = load.timescale()    
    
    intrack_offset_array = np.linspace(-10,10,20)
    offtrack_offset_array = np.linspace(-0.5,0.5, 20)

    all_combinations = [(a,b) for a in intrack_offset_array for b in offtrack_offset_array]

    waterfall = np.zeros((20, 20, args.channels))

    for f in tqdm(range(args.channels)):
         if args.debug:
             print("working on channel {}".format(f))
         global utc_array, cube
         utc_array, cube = getCube(f)

         with multiprocessing.Pool(processes=20) as pool:
             results = pool.starmap(worker, all_combinations)

         results = np.array(results)
         waterfall[:,:,f] = results.reshape((20,20))




    #for f in range(args.channels):
    #    if args.debug:
    #        print("working on channel {0}".format(f))
    #    utc, cube = getCube(f)
    #  
    #    ## mpi rest of the work
    #    with multiprocessing.Pool(processes=20) as pool:
    #        results = pool.starmap(worker, all_combinations)



        


if __name__ == "__main__":
 
    parser = ArgumentParser("satStack3D", description="stacks and integrates over over the predicted path of the satellite")
    parser.add_argument("--obs",required=True,type=int, help="the observation id")
    parser.add_argument("--noradid",required=True,type=int, help="the norad id of the satellite")
    parser.add_argument("--debug",type=bool, default=False, help="run script in debug mode?")
    parser.add_argument("--timeSteps",default=54,type=int, help="the number of timesteps")
    parser.add_argument("--channels",default=768,type=int, help="the number of fine channels")
    parser.add_argument("--user", required=True, help="the user name for spacetrack.org")
    parser.add_argument("--passwd",required=True, help="the password for spacetrack.org")
    parser.add_argument("--midName", default="2m", help="the mid name of the fits file")
    args = parser.parse_args()

    if args.debug:
        print("running script in debug mode")

    main()
