from scipy.optimize import fmin_cobyla
import sys, os, subprocess, numpy

#P, E, L, w = 1000.0, 69e9, 0.5, 0.1 # N, Pa, m, m

# all values are scaled by 1000 for the mesh to be of decent size

fileName = 'optimizeTest.txt'
resultName = '../TestsResults/opt.txt'
resultsName = "../TestsResults/cobylaSpringResults.txt"

results_write = open(resultsName, 'w')

def g0(x):
    if x > -.04:
      return x+.04
    return x

def objective(x):
    print 'g0 of ', x
    xpos = 0.0
    if type(x) is numpy.ndarray:
      xpos = x[0]
    else:
      xpos = x
    #if (height < 0):
    #  print 'Displacement for', height, 'is ::', "NaN", '\n' 
    #  return -1000000
    print 'Calculating for xpos:', xpos, '\n'
    
    file_write = open(fileName, 'w')
    file_write.write(fileName + ".scad "+str(xpos*1000)+" 20.0")
    file_write.close()
    subprocess.check_output(['python', 'pipeline.py', '--template', 'templateBSpring.py', '--batch', fileName, '--sConfig', 'slic3rConfig.ini', '-c'])
    # read results from file and return those
    opt = open(resultName)
    for line in opt.readlines():
      curLine = line.strip().split(' ')
      print 'Displacement for', xpos, 'is ::', float(curLine[0]) / 1000, '\n' 
      print 1e-4 - float(curLine[0]) / -1000
      results_write.write(str(xpos) + " " + str(float(curLine[0]) / -1000) + "\n")
      return (float(curLine[0]) / -1000)

    results_write.write(str(xpos) + " NEWTON ERROR\n")

    return 10000.0
    # Displacement constraint if no result specified
    #height = x
    #I = w * height**3 / 12 # m^4
    #tip_disp = (P * L**3)/(3*E*I)
    #print "Displacement", 1e-4 - tip_disp
    #results_write.write(str(height) + " " + str(tip_disp) + "\n")
    #return 1e-4 - tip_disp # max(disp) < 1e-3 m (1 mm)

def g1(x):
    # height > 0.01 m (10 mm)
    if x > .04:
      return -x
    return x+.04
    
def g2(x):
    # height < 0.5 m (500 mm)
    return 100
    
h0 = 0.0 # 20 mm
constraints = [g0, g1, g2]
h_opt = fmin_cobyla(objective, h0, constraints, rhoend=1e-6, maxfun=100, catol=1e-6)
results_write.write(str(h_opt))
#results_write.write(objective(h_opt))
results_write.close()
print h_opt, objective(h_opt), g0(h_opt)
