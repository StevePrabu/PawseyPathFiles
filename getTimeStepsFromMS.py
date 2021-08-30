#!/usr/bin/env python
from __future__ import division, print_function
from argparse import ArgumentParser
from casacore.tables import table, taql
import numpy as np

def listTimeStamps(ms):
    """
    function lists out all unique time-stamps in a ms set
    """
    mstime = ms.getcol("TIME")
    unique_times = np.unique(mstime)
    diff = unique_times[1:] - unique_times[:-1]
    output_times = [] 
    
    for t in unique_times:
        output_times.append(t)
        
    print("total of {} timesteps found in ms. integration time {}".format(len(output_times), np.mean(diff)))

    ## write info to file
    file1 = open("tmp.txt", "w")
    file1.write("TIMESTEPS={}\n".format(len(output_times)))
    file1.write("INTIME={}".format(np.mean(diff)))
    file1.close()

    return output_times


def main(args):
    ms = table(args.ms, readonly=True)
    time_steps = listTimeStamps(ms)
    ms.close()


if __name__ == "__main__":
    parser = ArgumentParser("timeSteps", description="finds number of timeSteps in ms")
    parser.add_argument("--ms", required=True, help="the measurement set name")
    args = parser.parse_args()

    main(args)
