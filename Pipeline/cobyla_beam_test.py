from scipy.optimize import fmin_cobyla
import sys, os, subprocess

P, E, L, w = 1000.0, 69e9, 0.5, 0.1 # N, Pa, m, m

fileName = 'optimize.txt'
resultName = 'optimizeResult.txt'

def objective(x):
    height = x #units in m
    volume = L * w * height
    return volume

def g0(x):
    height = x
    # write data to file to be used in pipeline
    file_write = open(fileName, 'w')
    file_write.write(fileName + ".scad 500 100 "+str(x))
    file_write.close()
    print subprocess.check_output(['python', 'pipeline', '--template', 'templateBeam.py', '--batch', fileName, '-c'])
    # read results from file and return those
    # to be implemented
    # Displacement constraint
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
