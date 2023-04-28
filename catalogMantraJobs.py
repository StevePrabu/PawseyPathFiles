from subprocess import call
import csv
from argparse import ArgumentParser

def getJobID(fileName):
    f = open(fileName)
    line = f.readline()

    while line:
        split_string = line.split(" ")
        if split_string[0] == "Submitted" and split_string[1] == "job:":
            jobid = int(split_string[2])
            break
        else:
            line = f.readline()

    return jobid


def main(args):

    jobid = getJobID(args.file)
    print("job id {} for obs {}".format(jobid, args.obs))

    with open("jobIdVsObsIdMapping.csv", "a") as vsc:
        thewriter = csv.writer(vsc)
        line = [jobid, args.obs]
        thewriter.writerow(line)

    

if __name__ == "__main__":
    parser = ArgumentParser("catalogMantraJobs", description="catalogs job ids and obs id for mantra-ray-client")
    parser.add_argument("--obs", required=True, type=int, help="the obs id")
    parser.add_argument("--file", required=True,help="the mantra-ray-client output file")
    args = parser.parse_args()

    main(args)