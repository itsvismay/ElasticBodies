from scipy.optimize import fmin_cobyla
import sys, os, subprocess, numpy

P, E, L, w = 1000.0, 69e9, 0.05, 0.01 # N, Pa, m, m

fileName = 'optimizeTest.txt'
resultName = '../TestsResults/opt.txt'

def objective(x):
    height = x
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

    # fail-safes
    if height <= 0.0:
      return -100
    
    file_write = open(fileName, 'w')
    file_write.write(fileName + ".scad 50 10 "+str(height*1000))
    file_write.close()
    subprocess.check_output(['python', 'pipeline.py', '--template', 'templateBeam.py', '--batch', fileName, '--sConfig', 'slic3rConfig.ini', '--ind', str(height* 1000), '-c'])
    # read results from file and return those
    opt = open(resultName)
    for line in opt.readlines():
      curLine = line.strip().split(' ')
      disp = float(curLine[0]) / 1000

      I = w * height**3 / 12
      tip_disp = (P * L**3)/(3*E*I)

      print 'Displacement for Height', height, 'is ::', disp
      print 'Analytical Disp for Height', height, 'is ::', tip_disp, '\n'

      return 1e-4 - (float(curLine[0]) / -1000)

    return -1000000    

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
