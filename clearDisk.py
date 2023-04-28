from subprocess import call
import csv

csv_files = ["hexEpoch1ObsSet1.csv", "hexEpoch1ObsSet2.csv"]

for fileName in csv_files:
    with open(fileName) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        for row in csv_reader:
            
            obs = int(row[0])
            print("working on obs {}".format(obs))
            syntax = "rm /astro/mwasci/sprabu/satellites/MWA-ORBITAL/processing/{}/*.npy".format(obs)
            bashExecute = call(syntax, shell=True)

            syntax2 = "rm /astro/mwasci/sprabu/satellites/MWA-ORBITAL/processing/{}/*6S*.fits".format(obs)
            bashExecute2 = call(syntax2, shell=True)