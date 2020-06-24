#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
from astropy.io import fits
from tqdm import tqdm
import csv
import os
from argparse import ArgumentParser

def main(obs):
     
    with open("missingFiles.csv", "w") as vsc:
        thewriter = csv.writer(vsc)
        thewriter.writerow(["t", "f"])
        
        for t in range(55):
            for f in range(768):
                
                if debug:
                    print("working on f {0} t {1}".format(f, t))
                if os.path.isfile(str(obs) + "-2m-" + str(t) + "-" + str(f).zfill(4)+ "-dirty.fits" ):
                    pass
                else:
                    if debug:
                        print("file not found " + str(obs) + "-2m-" + str(t) + "-" + str(f).zfill(4)+ "-dirty.fits")
                    thewriter.writerow([t, f])


             
            
            


if __name__ == "__main__":
    parser = ArgumentParser("findMissingFile", description="makes a list of missing fits images")
    parser.add_argument("--obs", required=True, help="the obs id")
    parser.add_argument("--debug", type=bool, default=False, help="runs script in debug mode")
    args = parser.parse_args()

    global debug
    debug = args.debug

    main(args.obs)
