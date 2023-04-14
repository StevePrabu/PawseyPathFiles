from __future__ import print_function, division
from argparse import ArgumentParser
from datetime import datetime
import csv

def main(args):

    ## read csv file and find max and min timestamps
    utc_array = []
    with open(args.config) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        for row in csv_reader:
            try:
                utc = datetime.strptime(row[4], '%Y-%m-%dT%H:%M:%S.%f')
            except:
                utc = datetime.strptime(row[4], '%Y-%m-%dT%H:%M:%S')
            utc_array.append(utc)

    ## fromat the datetime to required format
    start_utc = min(utc_array).strftime('%Y/%m/%d/%H:%M:%S.%f')
    end_utc = max(utc_array).strftime('%Y/%m/%d/%H:%M:%S.%f')

    print("time range to split {} - {}".format(start_utc, end_utc))

    ## casa split
    split(vis=args.inputMS, outputvis=args.outputMS,\
     datacolumn='all', timerange="{}~{}".format(start_utc, end_utc))

if __name__ == "__main__":
    parser = ArgumentParser("splitMS", description="splits the ms in time")
    parser.add_argument("--config", required=True, help="file that contains phase tracking info")
    parser.add_argument("--inputMS", required=True, help="the name of the input ms")
    parser.add_argument("--outputMS", required=True, help="the name of the output ms")
    args = parser.parse_args()
    main(args)