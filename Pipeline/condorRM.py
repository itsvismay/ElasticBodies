import sys, os, subprocess

val = int(sys.argv[1])
val2 = int(sys.argv[2])

for i in range(val, val2+1):
  print i
  try:
    subprocess.check_output(["condor_rm", str(i)])
  except:
    blah = 0
