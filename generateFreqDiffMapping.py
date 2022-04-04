#!/usr/bin/env python
from __future__ import print_function, division
from argparse import ArgumentParser
from casacore.tables import table, taql
import pandas as pd
import numpy as np

defaultPathFM = "/astro/mwasci/sprabu/path/PawseyPathFiles/FMinWA.txt"
"""
note all freq in this scipt is in MHz
"""


def getChannelsAndBandwidth(ms):
    """
    function finds the no. of channels, bandwidth and freq from ms
    """
    msfreq = ms.SPECTRAL_WINDOW.getcell("CHAN_FREQ",0)
    bandwidth = round(msfreq[1] - msfreq[0])/10**6
    noChannels = len(msfreq)
    print("the fine channel bandwidth is {}MHz and no. of channels {}".format(bandwidth, noChannels))
    
    return noChannels, bandwidth, msfreq

def indentifyFineChannels2Flag(noChannels, bandwidth, msfreq):
    
    edge_flag_bandwidth = 0.08
    centre_flag_bandwidth = 0.08
    
    channels_per_cc = int(noChannels/24) ## channels per coarse channel
    no_edge_channels2flag = int(edge_flag_bandwidth/bandwidth)
    no_cent_channels2flag = int(centre_flag_bandwidth/bandwidth)
    
    ## fine coarse channel index to flag
    flag_channels = []
    for i in range(no_edge_channels2flag):
        flag_channels.append(i)
        flag_channels.append((channels_per_cc-1)-i)
    flag_channels.append(int((channels_per_cc-1)/2))
    flag_channels.append(int((channels_per_cc-1)/2 +1))

    print("flag channel indexes {}".format(flag_channels))

    return flag_channels
    
def getFMfreq(args):
    """
    gets the fm freq
    """

    def getfreq(line):
        split_line = line.split(" ")
        for string in split_line:
            try:
                freq = float(string)
                if freq < 50:
                    continue
                break
            except:
                continue
        return freq
    
    def getPower(line, freq):
        split_line = line.split(" ")
        for ind in range(len(split_line)):
            string = split_line[ind]
            if string == freq:
                continue
            else:
                try:
                    if float(split_line[ind+1]) < 0:
                        power = split_line[ind]
                        break
                except:
                    continue
        try:
            return float(power)
        except:
            return float(power[:-1])*1000

    fm_freq_array = []
    f = open(args.path2transmitters)
    line = f.readline()
    
    while line:
        freq = getfreq(line)
        power = getPower(line, freq)
        if freq not in fm_freq_array and power >= args.powerCutoff:
            fm_freq_array.append(freq)
        
        line = f.readline()

    return fm_freq_array


def checkIfContainsFM(freq, bandwidth, fm_freq_array):
    fm_bandwidth = 0.2
    fm_half_band = fm_bandwidth/2
    half_bandwidth = bandwidth/2
    freq = freq/10**6

    for fm_freq in fm_freq_array:
        ## fine channel fully inside fm channel
        if freq+half_bandwidth <= fm_freq + fm_half_band \
        and freq-half_bandwidth >= fm_freq - fm_half_band:
            return True
        # fine channel partially inside fm channel
        if freq+half_bandwidth > fm_freq - fm_half_band \
        and freq-half_bandwidth < fm_freq - fm_half_band:
            return True
        if freq+half_bandwidth > fm_freq + fm_half_band \
        and freq-half_bandwidth < fm_freq + fm_half_band:
            return True
        
    return False


def main(args):

    if args.path2transmitters == defaultPathFM:
        print("external file for FM transmitters not provided")
        print("using default file found at {}".format(defaultPathFM))
    else:
        print("using fm transmitter file found at {}".format(args.path2transmitters))

    ## find no. of fine channels and bandwidth
    ms = table(args.ms, readonly=True)
    noChannels, bandwidth, msfreq = getChannelsAndBandwidth(ms)
    flag_channels = indentifyFineChannels2Flag(noChannels, bandwidth, msfreq)

    ## def get fm transmitters
    fm_freq_array = getFMfreq(args)

    print("found {} fm channels that satisfy conditions".format(len(fm_freq_array)))
    
    ## find channels not affected by FM
    not_affected_fine_channel = []
    not_affected_fine_channel_index = []
    MWA_freq_affected = []

    cf = 0
    for i in range(noChannels):
        if checkIfContainsFM(msfreq[i], bandwidth, fm_freq_array):
            MWA_freq_affected.append(True)
        else:
            MWA_freq_affected.append(False)
        
        if checkIfContainsFM(msfreq[i], bandwidth, fm_freq_array) ==  False and cf not in flag_channels:
            not_affected_fine_channel.append(msfreq[i])
            not_affected_fine_channel_index.append(i)
                        
        cf += 1
        if cf == int(noChannels/24):
            cf = 0

    print("{}/{} channels not affected by fm".format(len(not_affected_fine_channel), noChannels))

    ## find nearest fine channel to diff
    diff_array = []
    diff_channel_index = []
    for i in range(noChannels):
        ref_array = np.array(not_affected_fine_channel_index).T - i
        abs_ref_array = np.abs(ref_array)
        abs_ref_array[abs_ref_array == 0] = 10000
        
        solution = np.where(abs_ref_array == np.min(abs_ref_array))[0]
        if len(solution) == 1:
            sol_arg = int(solution)
        else:
            sol_arg = int(solution[0])
        diff_channel_index.append(not_affected_fine_channel_index[sol_arg])
        diff_array.append(np.abs(i - not_affected_fine_channel_index[sol_arg]))

    print("max channel diff is {} fine channels".format(np.nanmax(diff_array)))

    ## create dataframe and save to disk
    d = {"mwaChannelIndex" : range(noChannels), "diffChannelIndex": diff_channel_index, "diff" : diff_array, "FM affected": MWA_freq_affected}
    df = pd.DataFrame(data=d)
    print(df)
    df.to_pickle("freqDiffMap.plk")
    
    

if __name__ == "__main__":
    parser = ArgumentParser("generateFreqDiffMapping", description="creates a map of frequencies to perform freq diffferencing")
    parser.add_argument("--ms", required=True, help="the measurement set")
    parser.add_argument("--path2transmitters", default=defaultPathFM, help="the txt file that contains the FM transmitters")
    parser.add_argument("--powerCutoff", type=float, default=10000, help="the power cutoff for fm transmitters (default=10k)")
    args = parser.parse_args()

    main(args)



