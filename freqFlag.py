from __future__ import print_function, division
import numpy as np
from argparse import ArgumentParser
from tqdm import tqdm
import matplotlib.pyplot as plt
from casacore import tables
from argparse import ArgumentParser


def getFlags(freq_data, args):
    rms = np.sqrt(np.nanmean(freq_data**2))
    plt.plot(freq_data)
    plt.hlines(rms*args.threshold, 0, 768)
    plt.savefig("flagFreqInfo.png")
    flag_pos = np.where(freq_data >= args.threshold*rms)
    return flag_pos




def main(arg):
    
    ## get data
    measurementSet = tables.table(str(args.obs) + ".ms", readonly=False)
    data = measurementSet.getcol(arg.col)
    Flagdata = measurementSet.getcol("FLAG")
    inputShape = data.shape
    freq_data = np.abs(np.nansum(data, axis=(0, 2)))
    
    flag_pos = getFlags(freq_data, args)
    
    ## do flag
    flag = np.empty((inputShape[0], inputShape[-1]), dtype=bool)
    flag.fill(1)

    if flag_pos[0].shape == (0,):
        print("not freq to flag")
    else:
        print("found freq to flag")
        for f in flag_pos[0]:
            Flagdata[:,f,:] = flag
        print("flagging")
        measurementSet.putcol("FLAG", Flagdata)
        measurementSet.flush()
        measurementSet.close()
        print("done")





if __name__ == "__main__":
    parser = ArgumentParser("freqFlagger", description="flaggs frequency channels")
    parser.add_argument("--obs", required=True, type=int, help="the observation id")
    parser.add_argument("--col", required=True, help="the col in ms to flag")
    parser.add_argument("--threshold", default=5, type=int, help="the sigma threshold to flag")
    args = parser.parse_args()

    main(args)
