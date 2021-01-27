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
from os import path
import json
from subprocess import call
import csv


def obtainTLE(noradid,refUTC):
    time1 = refUTC + timedelta(hours=-24*2)
    time2 = refUTC 
    day1, month1, year1 = str(time1.day).zfill(2), str(time1.month).zfill(2), str(time1.year)
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



def getCutOff(distance):
    ## using nea-field equation from the xiang's metoer paper
    wavelength = 3.059106714 ## lambda at 98MHz
    d = np.sqrt(wavelength*distance/2)
    return d


def main(args):
  
    ## load metafits
    hdu = fits.open(args.metafits) 
    duration = hdu[0].header["EXPOSURE"]
    quack = float(hdu[0].header["QUACKTIM"])
    try:
        start_utc = datetime.strptime(hdu[0].header["DATE-OBS"], '%Y-%m-%dT%H:%M:%S.%f')
    except:
        start_utc = datetime.strptime(hdu[0].header["DATE-OBS"], '%Y-%m-%dT%H:%M:%S')
    start_utc =  start_utc + timedelta(seconds=quack) ## update start time for quack time
    pointing_ra = float(hdu[0].header["RA"])
    pointing_dec = float(hdu[0].header["DEC"])


    ## get tle 
    line1, line2, line3 = obtainTLE(args.noradid, start_utc)
    print(line2)
    print(line3)
    ts = load.timescale(builtin=True)
    mwa = Topos("26.701276778 S", "116.670846137 E", elevation_m= 377.827) 
    sat = EarthSatellite(line2, line3, line1, ts)

    ## create delta array (in seconds)
    search_timeSteps, search_ra, search_dec = [], [], [] ## contains the timesteps with satellite signal
    baseline_cutoff = []
    time_array = []
    for timeStep in range(int(duration/args.integration)):

        ### skip the first 2 timeSteps and the last timeStep
        if timeStep in [0, 1, int(duration/args.integration)-1]:
            continue
   
        local_utc = start_utc + timedelta(seconds=timeStep*args.integration)

        local_ts = ts.utc(local_utc.year, local_utc.month, local_utc.day, local_utc.hour, local_utc.minute, local_utc.second + local_utc.microsecond/1000000.0)

        sat.at(local_ts)
        difference = sat - mwa
        topocentric = difference.at(local_ts)
        ra, dec, distance = topocentric.radec()
        ra_deg, dec_deg = np.degrees(ra.radians), np.degrees(dec.radians)
        dist = np.sqrt((ra_deg- pointing_ra)**2 + (dec_deg - pointing_dec)**2)
        cutOff = getCutOff(distance.m)
        print("ra {} dec {} timeStep {} dist {} los dist {} cutoff {}".format(ra, dec, timeStep, dist, distance.m, cutOff))

        if dist < args.searchRadius:
            ## update ra and dec in chgcentre format
            ra = str(ra).replace(" ", "")
            dec = str(dec).replace(" ", "").replace('"', 's').replace("'", "m").replace("deg", "d")
            search_timeSteps.append(timeStep)
            search_ra.append(ra)
            search_dec.append(dec)
            baseline_cutoff.append(cutOff)
            time_array.append(str(local_utc))



    if debug:
        print("searching for sat in timeSteps {}".format(search_timeSteps))

    with open(str(args.obs)+ "-" + str(args.noradid) + ".csv", "w") as vsc:
        thewriter = csv.writer(vsc)
        for t, ra, dec, c, ut in zip(search_timeSteps,search_ra, search_dec, baseline_cutoff, time_array):
            line = [t, ra, dec, c, ut, "dummy"]
            thewriter.writerow(line)
   
    ## do the imaging
    #for t, ra, dec in zip(search_timeSteps, search_ra, search_dec):
    #    
    #    if debug:
    #        print("working on timeStep {}".format(t))

    #    ## change centre
        
    #    chg_synt = "singularity exec /pawsey/mwa/singularity/wsclean/wsclean_2.9.2-build-2.sif chgcentre {}.ms {} {}".format(args.obs, ra, dec)

    #    bashExecute = call(chg_synt, shell=True)
   
    #    ## image
    #    wsclean_synt = "singularity exec /pawsey/mwa/singularity/wsclean/wsclean_2.9.2-build-2.sif wsclean -name {0}-2m-{1} -size 100 100 -scale 5amin -interval {1} {2} -channels-out 768 -weight natural -abs-mem 10 {0}.ms".format(args.obs, t, t+1)

    #    bashExecute = call(wsclean_synt, shell=True)


        

        
        
        

    
    



if __name__ == "__main__":
    parser = ArgumentParser("track", description="finds timesteps where satellite is within FOV")
    parser.add_argument("--obs",required=True, type=int,help="the obs id")
    parser.add_argument("--metafits",required=True, help="the observation metafits file")
    parser.add_argument("--noradid",required=True, help="the noradid of the satellite")
    parser.add_argument("--integration",default=2,type=int,help="the image integration time")
    parser.add_argument("--searchRadius",default=18,type=float,help="the distance from pointing centre to search for satellite")
    parser.add_argument("--user",required=True, help="the user name for spacetrack.org")
    parser.add_argument("--passwd",required=True, help="the password for spacetrack.org")
    parser.add_argument("--debug",default=False,type=bool, help="run script in debug mode?Default=False")
    args = parser.parse_args()

    global debug, query
    debug = args.debug
    query = st.SpaceTrackClient(args.user, args.passwd)

    main(args)



