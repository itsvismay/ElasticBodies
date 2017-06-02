import sys, os, subprocess

#name = "Individual_"
name = "TwoVarCurve_"
ind = 0
done = 0

completedInd = []

for ind in range(0, 500):
  try:
    f = open(name+str(ind)+'/curvePrepedRemesh/position.txt')
    lines = f.readlines()
    #print ind, len(lines)
    #if len(lines) == 1602:
    if len(lines) >= 750:
      done = done + 1
      print ind
      completedInd.append(ind)
    #else:
    #  print ind
  except:
    blah = 0
    #print 'Not Found'

#print "Done:", done

#print completedInd

#subprocess.check_output(["ls", "/u/zmisso/"])

for ind in completedInd:
  print ind
  subprocess.check_output(["cp", name+str(ind)+"/curvePrepedRemesh/position.txt", "/u/zmisso/finished/"+name+str(ind)+"_Position.txt"])
  subprocess.check_output(["cp", name+str(ind)+"/points.txt", "/u/zmisso/finished/"+name+str(ind)+"_Points.txt"])
