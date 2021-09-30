#!/usr/bin/env python
from __future__ import print_function, division
import csv
from argparse import ArgumentParser
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def readDetectionFile(fileName):
    """
    extracts the x, xerr, y, yerr, and time information from csv file
    """
    x_array, y_array = [], []
    t_array = []
    with open(fileName) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")

        for row in csv_reader:
            x,x_err,y,y_err,RA,RA_err,DEC,DEC_err,UTC,timeStep = row

            if x == "x":
                continue
            else:

                x_array.append(float(x))
                y_array.append(float(y))
                
                t_array.append(float(timeStep))

    return x_array, y_array, t_array


def main(args):

    ## read detection file
    x_array, y_array, t_array = readDetectionFile(args.file)
    x_array = np.array(x_array)
    y_array = np.array(y_array)
    t_array = np.array(t_array)


    fig, ax = plt.subplots(figsize=(10, 10))
    plt.subplot(221)
    plt.scatter(t_array, x_array)
    plt.xlabel("ref time")
    plt.ylabel("x")

    plt.subplot(222)
    plt.scatter(t_array, y_array)
    plt.xlabel("ref time")
    plt.ylabel("y")

    dx_array = (x_array[1:] - x_array[:-1])/(t_array[1:] - t_array[:-1])
    dy_array = (y_array[1:] - y_array[:-1])/(t_array[1:] - t_array[:-1])

    plt.subplot(223)
    plt.scatter(t_array[1:], dx_array)
    plt.xlabel("ref time")
    plt.ylabel("dx/dt")

    plt.subplot(224)
    plt.scatter(t_array[1:], dy_array)
    plt.xlabel("ref time")
    plt.ylabel("dy/dt")    

    plt.savefig("validationAngularMeasurements{}.png".format(args.norad))





if __name__ == "__main__":
    parser = ArgumentParser("validate angular measurements", description="validates the extracted angular measurements")
    parser.add_argument("--file",required=True, help="the angular measurements of the pass")
    parser.add_argument("--norad",required=True, type=int,help="the noradid")
    args = parser.parse_args()

    main(args)