from scipy.optimize import fmin_cobyla
import sys, os, subprocess

P, E, L, w = 1000.0, 69e9, 0.5, 0.1 # N, Pa, m, m

# all values are scaled by 100 for the mesh to be of decent size

fileName = 'optimizeTest.txt'
resultName = '../TestResults/opt.txt'

def objective(x):
    height = x / 100 #units in m
    volume = L * w * height
    return volume

def g0(x):
    height = x
    print 'Calculating for Height:', height, '\n'
    # write data to file to be used in pipeline
    file_write = open(fileName, 'w')
    file_write.write(fileName + ".scad 50 10 "+str(x*100))
    file_write.close()
    print subprocess.check_output(['python', 'pipeline.py', '--template', 'templateBeam.py', '--batch', fileName, '--sConfig', 'slic3rConfig.ini', '-c'])
    # read results from file and return those
    opt = open(resultName)
    for line in opt.readlines():
      curLine = line.strip().split(' ')
      print 'Displacement for', height, 'is ::', float(curLine[0]) / 100, '\n' 
      return float(curLine[0]) / 100

    
    # Displacement constraint if no result specified
    #I = w * height**3 / 12 # m^4
    #tip_disp = (P * L**3)/(3*E*I)
    #return 1e-4 - tip_disp # max(disp) < 1e-4 m (0.1 mm)

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
