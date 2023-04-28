from subprocess import call
import csv
import subprocess
import datetime
import time


with open("obs2del.csv") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=",")
    for row in csv_reader:
        obs = int(row[0])
        print("working on obs {}".format(obs))
        syntax = "rm -r /astro/mwasci/sprabu/satellites/MWA-ORBITAL/processing/{}".format(obs)
        bashExecute = call(syntax, shell=True)
        

