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
import matplotlib.pyplot as plt
import os.path
from os import path
import json
from warnings import filterwarnings
filterwarnings("ignore")
import sys

def eprint(*args, **kwargs):
    """
    prints function in stderr so that the print is logged
    """
    print(*args, file=sys.stderr, **kwargs)


def getTLE(t):
    """
    gets the catalog of objects for the requested date range

    Parameters
    ----------
    t       : timestep for reference date (can be any)

    Returns
    -------
    catalog : the catalog of obtained objects
    entries : number of objects in catalog
    """

    time1 = startUTC + timedelta(hours=-24)
    time2 = time1 + timedelta(hours=24) ## obtain tle for objects updated within two weeks
    day1, month1, year1 = str(time1.day).zfill(2), str(time1.month).zfill(2), str(time1.year)
    day2, month2, year2 = str(time2.day).zfill(2), str(time2.month).zfill(2), str(time2.year)
    custom_name = year1 + "-" + month1 + "-" + day1 + "__" +  year2 + "-" + month2 + "-" + day2
    entries = 0

    if path.exists("../../catalogs/TLE_catalog" + custom_name + ".txt") == False:
                
        if debug:
            eprint("requesting file from server")

        date_range = ops.make_range_string(year1 + "-" + month1 + "-" + day1, year2 + "-" + month2 + "-" + day2)
        #result = query.tle_query(epoch=date_range)
        result = query.tle_publish_query(publish_epoch=date_range)

        ## write catalog to file
        with open("../../catalogs/TLE_catalog" + custom_name + ".txt", "w") as outfile:
            json.dump(result.json(), outfile)
        
        entries = len(result.json()[:])
        catalog = result.json()

    else:

        if debug:
            eprint("tle file found on disk. Not downloading.")

        with open("../../catalogs/TLE_catalog" + custom_name + ".txt") as json_file:
            result = json.load(json_file)

        entries = len(result[:])
        catalog = result
            
    return catalog, entries


def getSatXY(line1, line2, line3, UT, mwa, ts):
    """
    determines of the satellite is within fov and if so returns t
    he x y pixel location of the satellite

    Paramters
    ---------
    line1   : line1 of tle
    line2   : line2 of tle
    line3   : line3 of tle
    UT      : UTC time
    mwa     : skyfield observer object
    ts      : skyfield time object

    Returns
    -------
    px                  : x location of satellite
    py                  : y location of satellite
    number_of_pixels    : search cone radius
    visible             : visible? (bool)
    """

    visible = False
    number_of_pixels = 0
    sat = EarthSatellite(line2, line3, line1, ts)

    # propogate satellite to utc
    time = ts.utc(UT.year, UT.month,  UT.day, UT.hour, UT.minute, UT.second)
    sat.at(time)
    difference = sat - mwa
    topocentric = difference.at(time)

    # determine angular and pixel space location of satellite
    ra, dec, distance = topocentric.radec()
    ra, dec = np.degrees(ra.radians), np.degrees(dec.radians)
    px, py = wcs.all_world2pix([[ra, dec]], 1)[0]

    ## check if satellite within image

    if 0 < px < imgSize and 0 < py < imgSize:

        if np.sqrt((px - (imgSize/2))**2 + (py -(imgSize/2))**2) < imgSize/3:
            LOS_range = distance.m

            if LOS_range < 2000000:
                visible = True

            radius = np.degrees(70000/LOS_range) # 25 km offset cone in orbit (can be made smaller) s=rxtheta
            number_of_pixels = radius/pixel_scale 

    return px, py, number_of_pixels, visible


def checkDetected(px,py,radius,data):
    """
    creates a mask of given radius and checks if we a 
    event match
    
    Paramters
    ---------
    px      : the x pixel location
    py      : the y pixel location
    radius  : the radius of search cone in pixels
    data    : the data from rfiseeker

    Returns
    -------
    detected    : detected?? (bool)
    """
    x = np.linspace(0, (imgSize-1), imgSize)
    y = np.linspace(0, (imgSize-1), imgSize)
    x, y = np.meshgrid(x, y)

    array_radius = np.sqrt((x - px)**2 + (y - py)**2)
    array_radius[array_radius > radius] = -10
    array_radius[array_radius != -10] = 1
    array_radius[array_radius != 1] = 0

    masked_data = array_radius * data
    
    if np.sum(masked_data) > 0:
        detected = True

    else:
        detected = False

    return detected

def filterDetections(detectionSummary):
    """
    filters events that do not fit the below critera
    -appear more than three times
    -have atleast 2 consec envents
    -no detect timespacing less than 3

    Paramters
    ---------
    detectionSummary   : the dictionary of the detected events

    Returns
    -------
    filteredSummary    : the output detection summary
    """



    ## remove objects will less than 4 occurances and has no consec detections
    filteredSummary = dict()
    for obj in detectionSummary:

        noradid = detectionSummary[obj]["norad"]
        total = detectionSummary[obj]["total"]
        timeSteps = detectionSummary[obj]["timeSteps"]

        if float(total) > 4:

            timeStep_array = timeSteps.split("-")

            for i in timeStep_array:

                if str(int(i) + 1) in timeSteps:
                    
                    filteredSummary[obj] = dict(timeSteps=timeSteps, total=total, norad=noradid)

    

    ## removes detection with more that 3 timestep diff
    all_timeSteps = []
    for obj in filteredSummary:

        total = filteredSummary[obj]["total"]
        timeSteps = filteredSummary[obj]["timeSteps"]

        output_timeSteps = []
        output_total = 0
        timeStep_array = timeSteps.split("-")

        for i in timeStep_array:

            if str(int(i)+1) in timeStep_array or str(int(i)-1) in timeStep_array or \
                str(int(i)+2) in timeStep_array or str(int(i)-2) in timeStep_array or \
                    str(int(i)+3) in timeStep_array or str(int(i)-3) in timeStep_array:
                
                output_timeSteps.append(str(i))
                all_timeSteps.append(int(i))
                output_total += 1

        filteredSummary[obj]["total"] = str(output_total)

        filteredSummary[obj]["timeSteps"] = "-".join(output_timeSteps)

    imaging_timeSteps = []

    ## sort out the timeSteps for imaging
    for i in all_timeSteps:

        imaging_timeSteps.append(str(i))
        imaging_timeSteps.append(str(i+1))

    return filteredSummary, set(imaging_timeSteps)

# def updateTLE_2_obtain_latest(catalog, noObjects):
#     """
#     goes throug the catalog and outputs only the latest events
#     """


def main(obs, t1, t2):
    
    timeSteps = np.arange(t1, t2+1)

    if debug:
        eprint("the selected timeSteps are " + str(timeSteps))

    ## get tle catalog
    catalog, noObjects = getTLE(timeSteps[0])
    
    if debug:
        eprint("obtained TLE for {0} objects".format(noObjects))

    # begin search
    ts = load.timescale()
    mwa = Topos("26.703319405555554 S", "116.91558083333334 E", elevation_m= 377.827)
    globalData = np.zeros((imgSize, imgSize))

    ## the below was used for plotting...can be removed
    sat_x = []
    sat_y = []

    detectionSummary = dict()    

    for t in timeSteps:

        if debug:
            eprint("searching on timeStep {0}".format(t))

        ## get and time dati for the timestep 
        try:
            hdu = fits.open("6Sigma1Floodfilllr14SigmaRFIBinaryMap-t" + str(t).zfill(4) + ".fits")
            hduSeed = fits.open("6Sigma1Floodfilllr14SigmaRFIBinaryMapSeed-t" + str(t).zfill(4) + ".fits")
            dataSeed = hduSeed[0].data
        except:
            eprint("input file for timeStep " + str(t) + " not found.\nSkipping timestep.")
            continue
        data = hdu[0].data
        time = datetime.strptime(hdu[0].header['DATE-OBS'][:-2], '%Y-%m-%dT%H:%M:%S')
        globalData += data

        ## plot
        ax = plt.subplot(1,1,1, projection=wcs)
        tempData = np.zeros((imgSize, imgSize))
        tempData = globalData
        tempData = np.ma.masked_where(tempData == 0 , tempData)
        cmap= plt.cm.Purples
        cmap.set_bad(color="white")
        plt.imshow(tempData, origin="lower", vmax=1, vmin=-1, cmap= cmap)
        plt.title("UTC " + str(time))
        plt.xlabel("RA (Degrees)")
        plt.ylabel("DEC (Degrees)")
        plt.grid()
        sats_searched = [] # saves all the objects searched. Helps avoid duplicates
        
        for satNo in tqdm(range(noObjects)):
            
            #line1 = catalog[satNo]["OBJECT_NAME"]
            line1 = "sat"
            line2 = catalog[satNo]["TLE_LINE1"]
            line3 = catalog[satNo]["TLE_LINE2"]
            norad = line2[2:7]

            if norad in sats_searched:
                continue

            else:

                ######## find the most recent entry ################
                # value_array = []
                # for satNo2 in range(noObjects):
                #     line2 = catalog[satNo2]["TLE_LINE1"]
                #     norad_temp = line2[2:7]
                #     if float(norad_temp) == float(norad):
                #         epoch = datetime.strptime(catalog[satNo2]["PUBLISH_EPOCH"],'%Y-%m-%d %H:%M:%S' )
                #         diff = abs((epoch-startUTC).total_seconds())
                #         value_array.append([diff, satNo2])
                        
                # value_array = np.array(value_array)
                # recent_arg = int(np.where(value_array[:,0] == min(value_array[:,0]) )[0])
                # line2 = catalog[recent_arg]["TLE_LINE1"]
                # line3 = catalog[recent_arg]["TLE_LINE2"]

                #####################################################

                sats_searched.append(norad)
                #line2 = catalog[satNo]["TLE_LINE1"]
                #line3 = catalog[satNo]["TLE_LINE2"]
                #norad = line2[2:7]


                x, y, radius, visible = getSatXY(line1, line2, line3, time, mwa, ts)

                if visible:
                    detected = checkDetected(x, y , radius, dataSeed)
                    sat_x.append(x)
                    sat_y.append(y)

                    if detected:
                        circle = plt.Circle((x, y), radius, fill=False, edgecolor="lime")

                        if norad in detectionSummary:

                            ## obtain  prev values
                            prev_timeSteps = detectionSummary[str(norad)]["timeSteps"]
                            prev_total = detectionSummary[str(norad)]["total"]
                            current_timeSteps = prev_timeSteps + "-" + str(t)
                            current_total = str(int(prev_total) + 1)

                            ## update values
                            detectionSummary[str(norad)]["timeSteps"] = current_timeSteps
                            detectionSummary[str(norad)]["total"] = current_total
                          

                        else:

                            entry = dict(timeSteps=str(t), total="1", norad=norad)
                            detectionSummary[str(norad)] = entry

                    else:
                        circle = plt.Circle((x, y), radius, fill=False, edgecolor="red")

                    ax.add_artist(circle)

        plt.scatter(sat_x, sat_y, marker="x", color="magenta",s=1)
        plt.savefig("img" + str(t).zfill(2) + ".png", dpi=300)
        plt.clf()

    if debug:

        eprint("all detected events\n" + str(detectionSummary))

    output_summary, imaging_timeSteps = filterDetections(detectionSummary)
    
    if debug:
        eprint("filtered events\n" + str(output_summary))
        eprint("imaging timeSteps\n" + str(imaging_timeSteps))
    
    ## write to file
    json.dump(output_summary, open("filtered_summary.json", "w"))
    json.dump(detectionSummary, open("detection_summary.json", "w"))
    f1 = open("t.txt", "w")
    for l in imaging_timeSteps:
        f1.write(str(l)+"\n")
    f1.close()
    
    #### "That's all Folks!"" ###        
    if debug:
        eprint(" That's all Folks! ")



            

if __name__ == "__main__":

    parser = ArgumentParser("satSearch", description="searches for satellites in given data")
    parser.add_argument("--obs", required=True, type=int,help="The obs id")
    parser.add_argument("--t1", required=True, type=int,help="the start timeStep")
    parser.add_argument("--t2", required=True, type=int,help="the last timestep")
    parser.add_argument("--user", required=True, help="User name for space-track.org")
    parser.add_argument("--passwd", required=True, help="Password for space-track.org")
    parser.add_argument("--debug",default=False, type=bool, help="run scirpt in debug mode")
    args = parser.parse_args()
    
    global debug
    debug = args.debug
    
    global query
    query = st.SpaceTrackClient(args.user, args.passwd)

    ## get header info and make them global
    hdu = fits.open("6Sigma1Floodfilllr14SigmaRFIBinaryMap-t" + str(args.t1).zfill(4) + ".fits")
    
    global wcs, imgSize, startUTC, pixel_scale
    wcs = WCS(hdu[0].header, naxis=2)
    imgSize = hdu[0].header["NAXIS1"]
    pixel_scale = hdu[0].header["CDELT2"]
    startUTC = datetime.strptime(hdu[0].header['DATE-OBS'][:-2], '%Y-%m-%dT%H:%M:%S')

    if debug:
        eprint("running in debug mode")

    main(args.obs, args.t1, args.t2)
