#!/usr/bin/env python
from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from skyfield.api import Topos, load, EarthSatellite
import spacetracktool as st
from spacetracktool import operations as ops
from os import path
import json
from tqdm import tqdm
from argparse import ArgumentParser

def getTLE(startUTC):
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
    debug = True
    startUTC = str2datetime(startUTC)
    query = st.SpaceTrackClient("steverajprabu@gmail.com", 'Qwertyuiop1234567890')
    time1 = startUTC + timedelta(hours=-24*2)
    time2 = time1 + timedelta(hours=24) ## obtain tle for objects updated within two weeks
    day1, month1, year1 = str(time1.day).zfill(2), str(time1.month).zfill(2), str(time1.year)
    day2, month2, year2 = str(time2.day).zfill(2), str(time2.month).zfill(2), str(time2.year)
    custom_name = year1 + "-" + month1 + "-" + day1 + "__" +  year2 + "-" + month2 + "-" + day2
    entries = 0

    if path.exists("TLE_catalog" + custom_name + ".txt") == False:
                
        if debug:
            print("requesting file from server")

        date_range = ops.make_range_string(year1 + "-" + month1 + "-" + day1, year2 + "-" + month2 + "-" + day2)
        #result = query.tle_query(epoch=date_range)
        result = query.tle_publish_query(publish_epoch=date_range)

        ## write catalog to file
        with open("TLE_catalog" + custom_name + ".txt", "w") as outfile:
            json.dump(result.json(), outfile)
        
        entries = len(result.json()[:])
        catalog = result.json()

    else:

        if debug:
            print("tle file found on disk. Not downloading.")

        with open("TLE_catalog" + custom_name + ".txt") as json_file:
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
    #time = ts.utc(UT.year, UT.month,  UT.day, UT.hour, UT.minute, UT.second)
    time = UT
    sat.at(time)
    difference = sat - mwa
    topocentric = difference.at(time)

    # determine angular and pixel space location of satellite
    ra, dec, distance = topocentric.radec()
    ra, dec = np.degrees(ra.radians), np.degrees(dec.radians)
    #print(ra)
    #px, py = wcs.all_world2pix([ra, dec]], 1)[0]
    pix_coords = wcs.all_world2pix(np.array([ra, dec]).T, 1)
    px, py = pix_coords.T

    ## check if satellite within image

    imgSize = 1400
    pixel_scale = 0.0833333333333333

    return px, py

def str2datetime(str_date):
    """
    converts datetime in string format to dateime object
    """

    try:
        UT = datetime.strptime(str(str_date), '%Y-%m-%d %H:%M:%S')

    except:
        UT = datetime.strptime(str(str_date), '%Y-%m-%d %H:%M:%S.%f')

    return UT



def main(args):
    
    global_data = np.zeros((args.imgSize,args.imgSize))
    foundWCS = False
    global wcs
    utc_array = []
    for i in range(args.t1, args.t2):
        i+= 1
        hdu = fits.open(str(args.prefix) + "SigmaRFIBinaryMapPeakFlux-t"+str(i).zfill(4)+".fits")
        data = hdu[0].data
        global_data+= data
        utc = datetime.strptime(hdu[0].header["DATE-OBS"],'%Y-%m-%dT%H:%M:%S.%f')
        utc_array.append(utc)
        if not foundWCS:
            wcs = WCS(hdu[0].header, naxis=2)
            foundWCS = True

    #print(utc_array)
    ax = plt.subplot(1,1,1, projection=wcs)
    value_limit  = 1
    plt.imshow(global_data, cmap=plt.cm.seismic,  origin="lower", vmin=-value_limit, vmax=value_limit)
    plt.colorbar()
    plt.grid()
    ts = load.timescale()
    mwa = Topos("26.703319405555554 S", "116.91558083333334 E", elevation_m= 377.827)
    ### do satellite plotting
    def createTsTimeVector(utc_array):
        year_array, month_array, day_array, hour_array, minute_array, second_array = [],[],[],[],[],[]
        for utc in utc_array:
            start_utc = str2datetime(utc)
            start_date, start_time = str(start_utc).split(" ")
            start_year, start_month, start_day = start_date.split("-")
            start_hour, start_minute, start_second = start_time.split(":")
            year_array.append(int(start_year))
            month_array.append(int(start_month))
            day_array.append(int(start_day))
            hour_array.append(int(start_hour))
            minute_array.append(int(start_minute))
            second_array.append(float(start_second))
        times = ts.utc(year_array, month_array, day_array, hour_array, minute_array, second_array)
        return times
    time_vector = createTsTimeVector(utc_array)
    catalog, noObjects = getTLE(utc_array[0])


    ### plot sat track
    sats_searched = [] 
    for satNo in tqdm(range(noObjects)):
        line1 = "sat"
        line2 = catalog[satNo]["TLE_LINE1"]
        line3 = catalog[satNo]["TLE_LINE2"]
        norad = line2[2:7]
        if norad in sats_searched:
            continue
        else:
            sats_searched.append(norad)
            x, y = getSatXY(line1, line2, line3, time_vector, mwa, ts)
            if np.any(np.isnan(x)) or np.any(np.isnan(y)):
                continue
            else:
                
                plt.plot(x[:-1],y[:-1], color="blue", alpha=0.1)
                
    plt.savefig(args.obs + "_combinedRFIPlot.png")






if __name__ == "__main__":
    parser = ArgumentParser("combinedRFIPlot", description="makes combined rfi plot")
    parser.add_argument("--t1", default=1, type=int, help="the start timestep")
    parser.add_argument("--t2", default=55, type=int, help="the last time step")
    parser.add_argument("--obs", required=True, help="the obs id")
    parser.add_argument("--hpc",default="pawsey", help="the name of the hpc used to process data")
    parser.add_argument("--file", default="fits", help="make plot using fits file or vo table??, default fits")
    parser.add_argument("--prefix", required=True, help="the file prefix to obtain the wcs info")
    parser.add_argument("--imgSize", default=1400, type=int, help="the image size")
    args = parser.parse_args()

    main(args)
