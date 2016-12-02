from scipy.optimize import fmin_cobyla
import sys, os, subprocess, numpy

P, E, l, h = 1000.0, 69e9, 0.5, 0.1 # N, Pa, m, m

# all values are scaled by 1000 for the mesh to be of decent size

fileName = 'optimizeTest.txt'
resultName = '../TestsResults/opt.txt'
resultsName = "../TestsResults/cobylaResultsWidth.txt"

results_write = open(resultsName, 'w')

def objective(x):
    width = x #units in m
    volume = width * l * h
    return volume

def g0(x):
    print 'g0 of ', x
    width = 0.0
    if type(x) is numpy.ndarray:
      width = x[0]
    else:
      width = x
    if (width < 0):
      print 'Displacement for', width, 'is ::', "NaN", '\n' 
      return -1000000
    print 'Calculating for Width:', width, '\n'
    
    #file_write = open(fileName, 'w')
    #file_write.write(fileName + ".scad 50 "+str(width*1000)+" 10")
    #file_write.close()
    #subprocess.check_output(['python', 'pipeline.py', '--template', 'templateBeam.py', '--batch', fileName, '--sConfig', 'slic3rConfig.ini', '-c'])
    # read results from file and return those
    #opt = open(resultName)
    #for line in opt.readlines():
    #  curLine = line.strip().split(' ')
    #  print 'Displacement for', width, 'is ::', float(curLine[0]) / 1000, '\n' 
    #  print 1e-4 - float(curLine[0]) / -1000
    #  results_write.write(str(width) + " " + str(float(curLine[0]) / -1000) + "\n")
    #  return 1e-4 - float(curLine[0]) / -1000

    
    # Displacement constraint if no result specified
    width = x
    I = width * h**3 / 12 # m^4
    tip_disp = (P * l**3)/(E*I*3)
    print "Displacement", 1e-4 - tip_disp
    results_write.write(str(width) + " FAILED " + str(tip_disp) + "\n")
    print "blah blah"
    return 1e-4 - tip_disp # max(disp) < 1e-3 m (1 mm)

def g1(x):
    # length > 0.01 m (10 mm)
    return x - 0.01
    
def g2(x):
    # length < 0.5 m (500 mm)
    return 0.5 - x
    
h0 = 0.02 # 20 mm
constraints = [g0, g1, g2]
h_opt = fmin_cobyla(objective, h0, constraints, rhoend=1e-6, maxfun=100, catol=1e-6)
results_write.write(str(h_opt))
#results_write.write(objective(h_opt))
results_write.close()
print h_opt, objective(h_opt), g0(h_opt)
