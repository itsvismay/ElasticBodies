from scipy.optimize import fmin_cobyla
import sys, os, subprocess, numpy

P, E, L, w = 1000.0, 69e9, 0.5, 0.1 # N, Pa, m, m

# all values are scaled by 100 for the mesh to be of decent size

fileName = 'optimizeTest.txt'
resultName = '../TestsResults/opt.txt'

def objective(x):
    height = x #units in m
    volume = L * w * height
    return volume

def g0(x):
    print 'g0 of ', x
    height = 0.0
    if type(x) is numpy.ndarray:
      height = x[0]
    else:
      height = x
    print 'Calculating for Height:', height, '\n'
    #print 'what', height[0]
    # write data to file to be used in pipeline
    #if type(height) is numpy.ndarray:
    #  print 'WHY PYSCI WHY'
    #  print height[0]
    #print type(height)
    
    file_write = open(fileName, 'w')
    file_write.write(fileName + ".scad 50 10 "+str(height*1000))
    file_write.close()
    subprocess.check_output(['python', 'pipeline.py', '--template', 'templateBeam.py', '--batch', fileName, '--sConfig', 'slic3rConfig.ini', '-c'])
    # read results from file and return those
    opt = open(resultName)
    for line in opt.readlines():
      curLine = line.strip().split(' ')
      print 'Displacement for', height, 'is ::', float(curLine[0]) / 1000, '\n' 
      print 1e-4 - float(curLine[0]) / -1000
      # the 1e-3 is wrong, i just tried to use it to offset the values being
      # returned by the simulation (they were larger than they should have been)
      return 1e-3 - float(curLine[0]) / -1000

    
    # Displacement constraint if no result specified
    I = w * height**3 / 12 # m^4
    tip_disp = (P * L**3)/(3*E*I)
    return 1e-4 - tip_disp # max(disp) < 1e-4 m (0.1 mm)

def g1(x):
    # height > 0.01 m (10 mm)
    return x - 0.01
    
def g2(x):
    # height < 0.5 m (500 mm)
    return 0.5 - x
    
h0 = 0.02 # 20 mm
constraints = [g0, g1, g2]
h_opt = fmin_cobyla(objective, h0, constraints, rhoend=1e-6, maxfun=100, catol=1e-6)
print h_opt, objective(h_opt), g0(h_opt)
