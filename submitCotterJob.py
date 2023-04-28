from subprocess import call
import csv
import subprocess
import datetime
import time


def submitJob():
    ## get num of running jobs
    result = subprocess.check_output("squeue -u sprabu -M garrawarla", shell=True)
    runningJobs = result.count("sprabu")

    ## get num of pending jobs
    result = subprocess.check_output("squeue -l -M garrawarla", shell=True)
    pendingJobs = result.count("PENDING")

    print("found {} jobs currently in queue and {} pending jobs. Time={}".format(runningJobs, pendingJobs, datetime.datetime.now()))

    if runningJobs < 40:
        return True
    elif pendingJobs == 0:
        return True

    else:
        return False


def getASVOJobID(obs):
    with open("jobIdVsObsIdMapping.csv") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        for row in csv_reader:
            getASVOJobID, lobs = row
            getASVOJobID, lobs = int(getASVOJobID), int(lobs)

            if lobs == obs:
                return getASVOJobID


def main():

    with open("hexSurveyEpoch4obs.csv") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        for row in csv_reader:
            obs = int(row[0])
            asvoJobID = getASVOJobID(obs)
            syntax = "./bin/obs_cotter.sh -o " + str(obs) + " -j " + str(asvoJobID) + " -c /astro/mwasci/sprabu/satellites/DUG-MWA-SSA-Pipeline/calibrationSols/FMSurvey4.bin"
            bashExecute = call(syntax, shell=True)
            print("submitting job for obs {}".format(obs))
            while not submitJob():
                time.sleep(5*60) ## sleep for 5mins

if __name__ == "__main__":
    main()