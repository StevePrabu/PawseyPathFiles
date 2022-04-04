#!/usr/bin/env python
from __future__ import print_function, division
from argparse import ArgumentParser


def main(args):

    f = open(args.file)
    line = f.readline()

    while line:
        split_string = line.split(" ")
        if split_string[0] == "Submitted" and split_string[1] == "job:":
            jobid = int(split_string[2])
            break
        else:
            line = f.readline()

    print(jobid)


if __name__ == "__main__":
    parser = ArgumentParser("get job id", description="gets mwa client job id")
    parser.add_argument("--file", required=True, help="logger file that contains output from mantra ray client")
    args = parser.parse_args()

    main(args)