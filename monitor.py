#!/usr/bin/env python
from __future__ import division, print_function
import psutil
from argparse import ArgumentParser
import time

def getUsage():
    cpu = psutil.cpu_percent()
    memory = psutil.virtual_memory().percent
    return cpu, memory, psutil.disk_io_counters().read_bytes, psutil.disk_io_counters().write_bytes

def main(args):

    the_file = open(args.name + ".txt", 'a')
    the_file.write("cpu,memory,read,write\n")
    while True:
        cpu, memory, read, write = getUsage()
        
        the_file.write("{},{},{},{}\n".format(cpu, memory,read,write))
        time.sleep(args.sleepTime*60)


if __name__ == "__main__":
    parser = ArgumentParser("monitor", description="monitors cpu and memory usage")
    parser.add_argument("--name", required=True, help="name of the file that contains the usage info")
    parser.add_argument("--sleepTime", default=1, type=float, help="the sleep time between every entry in minutes")
    args = parser.parse_args()

    main(args)
