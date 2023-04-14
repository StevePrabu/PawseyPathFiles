#!/usr/bin/env python
from __future__ import division, print_function
from argparse import ArgumentParser
import time
import os

"""
the script checks if a file has a read lock. This script was motivated from the multiple 
ms read lock issuse faced by the ShiftStack pipeline on SpaceFest data.
"""

def main(args):

    while not os.access(args.file):

        ## if the file cannot be acessed, wait 5mins and check again
        print('file cannot be acessed at {}. Sleeping 5mins...'.time.ctime())
        time.sleep(5*60)

    print('file can be accessed.')


if __name__ == "__main__":
    parser = ArgumentParser('checkReadLock', description="checks if a file has read lock")
    parser.add_argument('--file', required=True, help="the input file")

    main(args)

