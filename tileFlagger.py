from __future__ import print_function, division
from casacore.tables import table, taql
import numpy as np
from argparse import ArgumentParser
#from tqdm import tqdm
import matplotlib.pyplot as plt
from subprocess import call

def main(args):

    mwa_tile_matrix = np.zeros((128, 128))
    mwa_tile_matrix[:] = np.nan

    input_ms = str(args.obs) + ".ms"
    ms = table(input_ms)

    for i in range(128):
        for j in range(128):

            if i == j: # no auto-correlations
                continue 
            
            try:
                temp = taql("select " + str(args.data_col) + " from $ms \
                    where ANTENNA1=" + str(i) + " and ANTENNA2=" + str(j))

                mwa_tile_matrix[i, j] = abs(np.sum(temp.getcol(str(args.data_col))))
            except:
                continue

    plt.imshow(mwa_tile_matrix, origin="lower", vmax=5* np.nanmean(mwa_tile_matrix), vmin=0)
    plt.colorbar()
    plt.savefig("tile" + str(args.obs) + ".png")
    rows, cols = np.where(mwa_tile_matrix >= args.threshold * np.nanmean(mwa_tile_matrix))

    set_rows, set_cols = set(rows), set(cols)

    flagged_tiles = []

    for row in set_rows:
        instance = len(np.array(np.where(rows == row))[0][:])
        if instance >= 10 and (row in flagged_tiles) == False:
            flagged_tiles.append(row)

    for col in set_cols:
        instance = len(np.array(np.where(cols == col))[0][:])
        if instance >= 10 and (col in flagged_tiles) == False:
            flagged_tiles.append(row)

    #print("tiles to be flagged " + str(flagged_tiles))


    if str(flagged_tiles) == "[]":
        print("no tiles to be flagged")
    else:
        print("tiles found to be flagged")
        for tile in flagged_tiles:
            print("flagging tile " + str(tile))
            synt1 = "flagged="+ str(tile)
            bashExecute1 = call(synt1, shell=True)
            synt2 = "echo $flagged | xargs flagantennae "+str(args.obs) + ".ms"
            bashExecute2 = call(synt2, shell=True)





    


if __name__ == "__main__":
    parser = ArgumentParser("tile flagger", description="flags tiles")
    parser.add_argument("--obs", required=True, type=int, help="the observation id")
    parser.add_argument("--data_col", default="DATA", help="the measurement set data column (DATA, CORRECTED_DATA)")
    parser.add_argument("--threshold", default=6, type=int, help="the threshold to flag")
    args = parser.parse_args()
    main(args)