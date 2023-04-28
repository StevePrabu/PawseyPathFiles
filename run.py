from subprocess import call
import csv
import subprocess


with open("events.csv") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=" ")
    for row in csv_reader:
        obs, norad = row[0], row[1]
        syntax = "./bin/obs_lightcurve.sh -o {} -n {} -c /astro/mwasci/sprabu/satellites/MWA-ORBITAL/hexSurveyEpoch1TLECatalog.txt".format(obs, norad)
        bashExecute = call(syntax, shell=True)
    
