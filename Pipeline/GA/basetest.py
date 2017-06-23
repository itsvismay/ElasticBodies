import sys, getopt, os, subprocess, time

timesToRun = 60

for i in range(timesToRun):
    print subprocess.check_output(["python", "ga.py"])
