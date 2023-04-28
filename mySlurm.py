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

    if runningJobs < 20:
        return True
    elif pendingJobs == 0:
        return True

    else:
        return False

def main():

    ## read csv file

    job_array, norad_array, obs_array = [], [], []

    with open("SpaceFest2022BSearchList.csv") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        for row in csv_reader:
            obs = int(row[0])
            norad = int(row[1])
            syntax = "./bin/obs_shiftstack.sh -o " + str(obs) + " -n " + str(norad) + " -i 0.5"
            #bashExecute = call(syntax, shell=True)
            #bashExecute = subprocess.run(syntax, stdout=subprocess.PIPE)
            bashExecute = subprocess.check_output(syntax, shell=True)
            jobid = bashExecute.split()[-1]
            print("submitting job for obs {}".format(obs))

            job_array.append(jobid)
            norad_array.append(norad)
            obs_array.append(obs)

            while not submitJob():
                time.sleep(5*60) ## sleep for 5mins

    with open('joblog.csv', 'w') as vsc:
        thewriter = csv.writer(vsc)
        line = ['obs' , 'norad', 'jobid']
        thewriter.writerow(line)

        for obs, norad, jobid in zip(obs_array, norad_array, job_array):
            line = [obs, norad, jobid]
            thewriter.writerow(line)
        


if __name__ == "__main__":
    main()
