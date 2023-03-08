#!/usr/bin/env python
from __future__ import division, print_function
import numpy as np
from casacore.tables import table
from astropy.io import fits
from argparse import ArgumentParser
from datetime import datetime, timedelta
from astropy.coordinates import EarthLocation, AltAz, SkyCoord
import matplotlib.pyplot as plt
import astropy.units as u
import time as tm
from multiprocessing import Pool
import csv
from astropy.time import Time
import sys
from scipy.ndimage import  rotate


def msTime2UTC(time):
    """
    converts from the wierd ms time format to usable time
    Parameters
    ----------
    time    : the wierd ms time
    Returs
    ------
    time    : useful time
    """

    time_local = time/(24*3600)
    tmp = Time(time_local, format='mjd')

    return tmp.datetime

def getRelevantRows(mstime, msnoRows, timestep):
    """
    get the rows from ms that correspond to the considered timestep
    """

    start = tm.time()

    converted_time = msTime2UTC(mstime)

    relevant_rows = []

    for i in range(msnoRows):
        t = str(converted_time[i])
        if str(t) != str(timestep):
            continue
        else:
            relevant_rows.append(i)

    ## quit if no relevant rows found
    if len(relevant_rows) == 0:
        print('no relevant time-steps found. Quitting...')
        sys.exit(0)

    if verbose:
        elapsed = tm.time() - start
        print('found {}/{} relevant ms indexes in {}s'.format(len(relevant_rows), msnoRows, elapsed))

    return relevant_rows


def readConfig(fileName):
    """
    reads the shift-stack config file
    """
    t_array, ra_array, dec_array = [], [], []
    dist_array, utc_array = [], []
    with open(fileName) as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',')
        for row in csvreader:
            t,ra,dec,dist,utc,_ = row
            t_array.append(t)
            ra_array.append(ra)
            dec_array.append(dec)
            dist_array.append(dist)
            utc_array.append(utc)

    return t_array, ra_array, dec_array, dist_array, utc_array


def reprojUVW(uvw, ra, dec, ra0, dec0):
    u, v, w = uvw.T
    newuvw = np.empty_like(uvw)
    
    newuvw[:, 0] = (
        u * np.cos(ra - ra0)
        + v * np.sin(dec0) * np.sin(ra - ra0)
        - w * np.cos(dec0) * np.sin(ra - ra0)
    )
    newuvw[:, 1] = (
        -u * np.sin(dec) * np.sin(ra - ra0)
        + v * (np.sin(dec0) * np.sin(dec) * np.cos(ra - ra0) + np.cos(dec0) * np.cos(dec))
        + w * (np.sin(dec0) * np.cos(dec) - np.cos(dec0) * np.sin(dec) * np.cos(ra - ra0))
    )
    newuvw[:, 2] = (
        u * np.cos(dec) * np.sin(ra - ra0)
        + v * (np.cos(dec0) * np.sin(dec) - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0))
        + w * (np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0))
    )
    return newuvw

def phaseRotate(data, oldw, neww, wavelengths):
    offset = -2j * np.pi * (neww - oldw)
    phase = np.empty_like(data)
    
    for row in range(data.shape[0]):
        tmp = offset[row]/wavelengths
        for pol in range(4):
            phase[row, :, pol] = tmp
            
    return data*np.exp(phase)



def changePhaseCentre(data, uvw, ra0_rad, dec0_rad, ra1_rad, dec1_rad, wavelengths):
    """
    reprojects the data to the new phase centre
    """
    if verbose:
        start = tm.time()
    newuvw = reprojUVW(uvw, ra1_rad, dec1_rad, ra0_rad, dec0_rad)
    newdata = phaseRotate(data, uvw[:,2], newuvw[:,2], wavelengths)
    
    if verbose:
        elapsed = tm.time() - start
        print('changed phase centre in {}s'.format(elapsed))


    return newuvw, newdata

def getUVCellSize(pixelScale, imgSize):
    """
    note pixelScale is in radians
    and we assume a square image of size imgSize x imgSize
    """
    
    cellSize = 1/(pixelScale*imgSize)
    
    return cellSize

    
def MPIimage(channel):
    """
    makes images. note imscale is in amin
    """
    pixelScaleArcmin = 1
    imgSize = 1024

    if verbose:
        start = tm.time()

    pixeScale = np.radians(pixelScaleArcmin/60)

    uvgrid_xx = np.zeros((imgSize, imgSize), dtype=complex)
    uvgrid_yy = np.zeros((imgSize, imgSize), dtype=complex)
    uvgrid_count = np.zeros((imgSize, imgSize), dtype=int)
    uv_cellSize = getUVCellSize(pixeScale, imgSize)
    

    c = 0

    for entry in relevant_rows:
        ## filter out auto corr.
        ant1 = msant1[entry]
        ant2 = msant2[entry]
    
        if ant1 == ant2:
            c += 1
            continue
        
        ## check if flagged
        if msflagrows[entry] or np.all(msflags[entry] == True):
            c += 1
            continue


        u,v,w = lensedUVW[c]
        
        u_lamba, v_lamba = u/wavelengths[channel], v/wavelengths[channel]

        if msflags[entry, channel, 0] or msflags[entry, channel, 3]:
            c += 1
            continue
        
        ## convert wavelengths to pixel no
        u_pix = int(u_lamba/uv_cellSize + imgSize/2)
        v_pix = int(v_lamba/uv_cellSize + imgSize/2)
        

        if np.isnan(lensedData[c, channel, 0]) or np.isnan(lensedData[c, channel, 3]):
            c += 1
            continue

        try:    
            uvgrid_xx[u_pix, v_pix] += lensedData[c, channel, 0]
            uvgrid_yy[u_pix, v_pix] += lensedData[c, channel, 3]
            uvgrid_count[u_pix, v_pix] += 1
        except:
            print('u pix {} v pix {}'.format(u_pix, v_pix))
        c += 1

    ## normalise the uv cel values
    uvgrid_avg_xx = np.zeros((imgSize, imgSize), dtype=complex)
    uvgrid_avg_yy = np.zeros((imgSize, imgSize), dtype=complex)
    for row in range(imgSize):
        for col in range(imgSize):
            count = uvgrid_count[row, col]
            if count == 0:
                continue
            uvgrid_avg_xx[row, col] = uvgrid_xx[row, col]/count
            uvgrid_avg_yy[row, col] = uvgrid_yy[row, col]/count

    img_xx = np.abs(np.fft.ifftshift(np.fft.ifft2(uvgrid_avg_xx)))
    img_yy = np.abs(np.fft.ifftshift(np.fft.ifft2(uvgrid_avg_yy)))
    imgLarg = rotate(np.sqrt(img_xx**2 + img_yy**2), 90, order=5, reshape=5)
    img = imgLarg[512-128:512+128, 512-128:512+128]

    ## save to file
    hdun = fits.PrimaryHDU(img)
    hdun.writeto('img-t{}-f{}.fits'.format(ind, channel), overwrite=True)

    ## save png
    plt.imshow(img, origin='lower')
    plt.colorbar()
    plt.title('t{}-f{}'.format(ind, channel))
    plt.savefig('img-t{}-f{}.png'.format(ind, channel))
    plt.clf()
    
    return 

def getLatLonAlt(x, y, z):
    """
    convert geocentric cartesian coordinates to lat lon and alt
    Parameters
    ----------
    x   : x coordinate
    y   : y coordinate
    z   : z coordinate
    Returns
    -------
    lat : latitude
    lon : longitude
    alt : altitude
    """

    pos = EarthLocation.from_geocentric(x,y,z, unit="m").to_geodetic()

    return pos.lat.degree, pos.lon.degree, pos.height

def getAltAz(x, y, z):
    """
    computes the alt-az of phase centre from geometrical centre of array
    Parameters
    ----------
    x       : geocetric mean x position of the array
    y       : geocetric mean y position of the array
    z       : geocetric mean z position of the array
    Returns
    -------
    alt     : the alt (deg) of phase centre
    az      : the az (deg) of phase centre
    """

    lat, lon, height = getLatLonAlt(x, y, z)
    pos = EarthLocation(lon=lon*u.deg, lat=lat*u.deg, height=float(str(height)[:-1])*u.m)
    coord = SkyCoord(ms_ra_deg, ms_dec_deg, unit=(u.deg,u.deg))
    coord.time = timestamp + timedelta(hours=pos.lon.hourangle)
    coord = coord.transform_to(AltAz(obstime=timestamp, location=pos))

    return np.degrees(coord.alt.rad), np.degrees(coord.az.rad)



def getFocusXYZ(focus):
    """
    returns the x,y,z coordinates of focal point 
    Parameters
    ----------
    focus   : the focal distance of the array (km)
    Returns
    -------
    x       : x coordinate of the focus
    y       : y coordinate of the focus
    z       : z coordinate of the focus
    """

    ## calcualte the Geocentric mean coordinates of the array
    mwa_x, mwa_y, mwa_z = np.mean(msantloc[:,0]), np.mean(msantloc[:,1]), np.mean(msantloc[:,2])

    ## calcualte alt-az of focus from mean position of array
    mwa_alt, mwa_az = getAltAz(mwa_x, mwa_y, mwa_z)

    ## calcualte x,y,z using focus and alt-az angles
    z = focus*1e3*np.cos(np.radians(mwa_alt))*np.cos(np.radians(mwa_az))
    x = focus*1e3*np.cos(np.radians(mwa_alt))*np.sin(np.radians(mwa_az))
    y = focus*1e3*np.sin(np.radians(mwa_alt))

    return x, y, z


def getAntXYZ(antName):
    """
    gets antenna's x,y,z coordinates in local reference frame.
    
    Note that the mean North, South, and Height values from 
    metafits does not coincide with (0,0,0). Hence, the returned
    positions are explicitely offseted by the mean value to force
    the centre of the array to be at (0,0,0). This does NOT change
    the calculated delay as it is a linear offset in the arbirarily 
    defined cartesian reference frame and the coordinates of the focal
    point also assumses (0,0,0) to the centre of the array
    Parameters
    ----------
    antName : the name of the antenna
    Returns
    -------
    tile_x  : x coordinate of the tile in local ref frame
    tile_y  : y coordinate of the tile in local ref frame
    tile_z  : z coordinate of the tile in local ref frame
    """
    
    ## mean offset to apply
    z_offset = np.mean(hdu[1].data['North'])
    y_offset = np.mean(hdu[1].data['Height'])
    x_offset = -np.mean(hdu[1].data['East'])
    
    # index of tile in metafits file
    index = np.where(hdu[1].data['TileName'] == antName)[0][0]
    
    ## updated x,y,z coordinates of the tile
    tile_z = hdu[1].data['North'][index] - z_offset
    tile_y = hdu[1].data['Height'][index] - y_offset
    tile_x = -hdu[1].data['East'][index] - x_offset
    
    return tile_x, tile_y, tile_z


def getPhase(ant1, ant2, uvw):
    """
    calculates near-field phase/delay for a baseline
    Parameters
    ----------
    ant1    : antenna 1 of the baseline being considered
    ant2    : antenna 2 of the baseline being considered
    uvw     : far field uvw coordinates of the baseline
    Returns
    -------
    phi     : the phase correction to apply to vis
    new_w   : the calculated near-field delay (or w-term)
    """

    antName1 = msantname[ant1]
    antName2 = msantname[ant2]
    u, v, w = uvw

    ## antenna positions in local cartesian reference frame
    tile1_x, tile1_y, tile1_z = getAntXYZ(antName1)
    tile2_x, tile2_y, tile2_z = getAntXYZ(antName2)

    ## calculate distance from antennas to focal point
    r1 = np.sqrt((tile1_x + focus_x)**2 + (tile1_y - focus_y)**2\
              + (tile1_z - focus_z)**2)
    r2 = np.sqrt((tile2_x + focus_x)**2 + (tile2_y - focus_y)**2\
              + (tile2_z - focus_z)**2)

    ## calculate near-field delay
    new_w = r2 - r1
    phi = new_w - w

    return phi, new_w


def LEOLens(newdata, newuvw):
    """
    applies near-field correction
    """
    
    if verbose:
        start = tm.time()



    lensedData = np.copy(newdata)
    lensedUVW = np.copy(newuvw)

    c = 0
    for entry in relevant_rows:
        ant1 = msant1[entry]
        ant2 = msant2[entry]

        if ant1 != ant2:
            uvw = lensedUVW[c]
            phi, w_new = getPhase(ant1, ant2, uvw)

            ## updated old w with new w
            lensedUVW[c][2] = w_new

            ## phase rotate the visibilities
            phase = -2j*np.pi*phi/wavelengths
            for pol in range(4):
                tmp = np.copy(lensedData[c,:,pol])
                lensedData[c,:,pol] = tmp*np.exp(phase)

    
        c += 1

    if verbose:
        elapsed = tm.time() - start
        print('applied near-field delays in {}s'.format(elapsed))

    return lensedData, lensedUVW


def main(args):

    t_array, ra_array, dec_array, dist_array, utc_array = readConfig(args.config)

    ## read ms
    if verbose:
        start = tm.time()

    global msant1, msant2, msflags, msflagrows, wavelengths,\
     relevant_rows, lensedUVW, lensedData, msantname, msantloc,\
      ms_ra_deg, ms_dec_deg, timestamp, t

    ## read the measurement set
    ms = table(args.ms, readonly=True, ack=False)
    ms_ra_rad, ms_dec_rad = ms.FIELD.getcol('PHASE_DIR')[0,0]
    ms_ra_deg, ms_dec_deg = np.degrees(ms_ra_rad), np.degrees(ms_dec_rad)
    msdata = ms.getcol('CORRECTED_DATA')
    msuvw = ms.getcol('UVW')
    mstime = ms.getcol('TIME')
    msant1 = ms.getcol('ANTENNA1')
    msant2 = ms.getcol('ANTENNA2')
    msfreq = ms.SPECTRAL_WINDOW.getcell('CHAN_FREQ', 0)
    msantloc = ms.ANTENNA.getcol("POSITION")
    msantname = ms.ANTENNA.getcol('NAME')
    msflags = ms.getcol('FLAG')
    msflagrows = ms.getcol('FLAG_ROW')
    msnoRows = len(msant1)
    ms.close()
    wavelengths = 299792458/msfreq

    if verbose:
        elapsed = tm.time() - start
        print('data loaded from ms. time elapsed {}s'.format(elapsed))

    global focus_x, focus_y, focus_z, ind

    for ind, ra, dec, dist, utc in zip(t_array, ra_array, dec_array, dist_array, utc_array):
        
        ra1_rad = SkyCoord(ra, dec, frame='icrs', unit=(u.hourangle, u.deg)).ra.rad
        dec1_rad = SkyCoord(ra, dec, frame='icrs', unit=(u.hourangle, u.deg)).dec.rad
        print('working on timestep {} ra {} dec {} utc {}'.format(ind, ra, dec, utc))

        ## copy the current ra dec to ms_ra/dec values (as it is used by getFocus())
        ms_ra_deg = np.copy(np.degrees(ra1_rad))        
        ms_dec_deg = np.copy(np.degrees(dec1_rad))

        try:
            timestamp = datetime.strptime(utc, '%Y-%m-%dT%H:%M:%S.%f')
        except:
            timestamp = datetime.strptime(utc, '%Y-%m-%dT%H:%M:%S')

        #get data slice
        relevant_rows  = getRelevantRows(mstime, msnoRows, str(timestamp))
        data_slice = np.copy(msdata[relevant_rows, :,:])
        uvw_slice = np.copy(msuvw[relevant_rows])
        
        ## change phase centre
        newuvw, newdata = changePhaseCentre(data_slice, uvw_slice, ms_ra_rad, ms_dec_rad, ra1_rad, dec1_rad, wavelengths)

        ### calculate the x,y,z coordinate of the focal point in local terrestrial reference frame
        focus_x, focus_y, focus_z = getFocusXYZ(float(dist))
        
        ### apply near-field delay corrections
        lensedData, lensedUVW = LEOLens(newdata, newuvw)

        ## perform imaging
        start = tm.time()
        p = Pool(args.cores)
        p.map(MPIimage, range(768))
        elapsed = tm.time() - start
        print('fine channel images done. time elapsed {}s'.format(elapsed))


    


if __name__ == "__main__":
    parser = ArgumentParser('pyShiftStack', description='does shift-stacking in python')
    parser.add_argument('--ms', required=True, help='input measurement set')
    parser.add_argument('--config', required=True, help='config file')
    parser.add_argument('--metafits', required=True,\
                help='the metafits of the observation')
    parser.add_argument('--verbose', type=bool, default=False,\
                help='runs scipt in verbose mode (default=False)')
    parser.add_argument('--cores', type=int, default=12,\
                help='no of cpus to use. default 12')
    args = parser.parse_args()

    global verbose, hdu
    verbose = args.verbose
    hdu = fits.open(args.metafits)

    if verbose:
        print('running Shift-Stacking in verbose mode')

    main(args)
