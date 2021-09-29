#!/usr/bin/env python
from __future__ import division, print_function
from argparse import ArgumentParser
from casacore.tables import table, taql
import numpy as np
from tqdm import tqdm
from astropy.table import Table
import pandas as pd
import time as tm
import matplotlib.pyplot as plt


def getTilePosition(reqTile, tbl):
    
    allTileNames = tbl["TileName"]
    tileNorth = tbl["North"]
    tileEast = tbl["East"]

    for tile, north, east in zip(allTileNames, tileNorth, tileEast):
        if tile == reqTile:
            return north, east


def main(args):

    ## load metafits table
    tbl = Table.read(args.metafits)
    
    ## load measurement set
    ms = table(args.ms, readonly=True)

    if debug:
        print("loading data from measurement set...",end="")
        start = tm.time()
    
    mstime = ms.getcol("TIME")
    msant1 = ms.getcol("ANTENNA1")
    msant2 = ms.getcol("ANTENNA2")
    msflags = ms.getcol("FLAG")
    msflagrows = ms.getcol("FLAG_ROW")
    msantname = ms.ANTENNA.getcol("NAME")
    ms.close()

    if debug:
        elapsed = tm.time() - start
        print("done. time taken {}s".format(elapsed))

    ## search for available tiles
    available_antennas = []
    available_tileNames = []
    available_north = []
    available_east = []

    for i in tqdm(range(len(mstime))):
    
        ant1 = msant1[i]
        ant2 = msant1[i]
        if np.all(msflags[i] == True) :
            continue

        if ant1 not in available_antennas:
            available_antennas.append(ant1)
            available_tileNames.append(msantname[ant1])
            north, east = getTilePosition(msantname[ant1], tbl)
            available_north.append(north)
            available_east.append(east)

        if ant2 not in available_antennas:
            available_antennas.append(ant2)
            available_tileNames.append(msantname[ant2])
            north, east = getTilePosition(msantname[ant2], tbl)
            available_north.append(north)
            available_east.append(east)

        
    ## create data frame
    if debug:
        print("a total of {}/128 tiles not flagged".format(len(available_antennas)))

    d = {'Tile available': available_antennas, 'Tile Name': available_tileNames,\
        'North': available_north, 'East': available_east}

    df = pd.DataFrame(d)
    
    df.to_pickle("availableTilesForObs{}.pkl".format(args.obs))



if __name__ == "__main__":
    parser = ArgumentParser("getTilesNotFlagged", description="gets the tiles not flagged")
    parser.add_argument("--ms", required=True, help="the measurement set name")
    parser.add_argument("--metafits", required=True, help="the metafits file for the observation")
    parser.add_argument("--obs", required=True, type=int, help="the observation id")
    parser.add_argument("--debug", default=False, type=bool, help="runs script in debug mode")
    args = parser.parse_args()

    global debug
    debug = args.debug

    main(args)

