from scipy.optimize import fmin_cobyla
import sys, os, subprocess, numpy

# all values are scaled by 1000 for the mesh to be of decent size

fileName = 'optimizeTest.txt'
resultName = '../TestsResults/opt.txt'
resultsName = "../TestsResults/cobylaSpringResults.txt"

results_write = open(resultsName, 'w')

### VALUES
#
# x[0]      x1
# x[1]      x2
# x[2]      y2
# x[3]      x3
# x[4]      y3
# x[5]      x4
# x[6]      y4
# x[7]      x5
# x[8]      w2
# x[9]      w3
# x[10]     w4
# x[11]     thk
# x[12]     hei
#

# min-bounds
def g0(x):
    if x[0] < -.04: # x1
      return -1
    if x[1] < -.04: # x2
      return -1
    if x[2] < -.04: # y2
      return -1
    if x[3] < -.04: # x3
      return -1
    if x[4] < -.04: # y3
      return -1
    if x[5] < -.04: # x4
      return -1
    if x[6] < -.04: # y4
      return -1
    if x[7] < -.04: # x5
      return -1
    if x[8] < 0.001: # w2
      return -1
    if x[9] < 0.001: # w3
      return -1
    if x[10] < 0.001: # w4
      return -1
    if x[11] < 0.001: # thk
      return -1
    if x[12] < 0.001: # hei
      return -1
    return 1

def objective(x):
    print 'g0 of ', x
    x1 = 0.0
    x2 = 0.0
    y2 = 0.0
    x3 = 0.0
    y3 = 0.0
    x4 = 0.0
    y4 = 0.0
    x5 = 0.0
    w2 = 0.0
    w3 = 0.0
    w4 = 0.0
    thk = 0.0
    hei = 0.0
    if type(x) is numpy.ndarray:
      x1 = x[0]
      x2 = x[1]
      y2 = x[2]
      x3 = x[3]
      y3 = x[4]
      x4 = x[5]
      y4 = x[6]
      x5 = x[7]
      w2 = x[8]
      w3 = x[9]
      w4 = x[10]
      thk = x[11]
      hei = x[12]
    else:
      x1 = x[0]
      x2 = x[1]
      y2 = x[2]
      x3 = x[3]
      y3 = x[4]
      x4 = x[5]
      y4 = x[6]
      x5 = x[7]
      w2 = x[8]
      w3 = x[9]
      w4 = x[10]
      thk = x[11]
      hei = x[12]
    
    print 'Calculating for values:', '\n'
    print 'X1:', x1
    print 'X2:', x2
    print 'Y2:', y2
    print 'X3:', x3
    print 'Y3:', y3
    print 'X4:', x4
    print 'Y4:', y4
    print 'X5:', x5
    print 'W2:', w2
    print 'W3:', w3
    print 'W4:', w4
    print 'THK:', thk
    print 'HEI:', hei
    
    file_write = open(fileName, 'w')
    file_write.write(fileName + ".scad "+str(x1*1000)+" "+str(x2*1000)+" "+str(y2*1000)+" "+str(x3*1000)+" "+str(y3*1000)+" "+str(x4*1000)+" "+str(y4*1000)+" "+str(x5*1000)+" "+str(w2*1000)+" "+str(w3*1000)+" "+str(w4*1000)+" "+str(thk*1000)+" "+str(hei*1000))
    file_write.close()
    subprocess.check_output(['python', 'pipeline.py', '--template', 'templateBSpring.py', '--batch', fileName, '--sConfig', 'slic3rConfig.ini', '-c'])
    # read results from file and return those
    opt = open(resultName)
    for line in opt.readlines():
      curLine = line.strip().split(' ')
      # not sure what this should be returning and testing for
      return (float(curLine[0]) / -1000)

    #results_write.write(str(xpos) + " NEWTON ERROR\n")

    return -10000.0

# max-bounds
def g1(x):
    # height > 0.01 m (10 mm)
    if x[0] > 0.04: # x1
      return -1
    if x[1] > 0.04: # x2
      return -1
    if x[2] > 0.04: # y2
      return -1
    if x[3] > 0.04: # x3
      return -1
    if x[4] > 0.04: # y3
      return -1
    if x[5] > 0.04: # x4
      return -1
    if x[6] > 0.04: # y4
      return -1
    if x[7] > 0.04: # x5
      return -1
    if x[8] > 0.01: # w2
      return -1
    if x[9] > 0.01: # w3
      return -1
    if x[10] > 0.01: # w4
      return -1
    if x[11] > 0.01: # thk
      return -1
    if x[12] > 0.1: # hei
      return -1
    return 1

# overlap    
def g2(x):
    if x[2] <= 0.0: # y2 <= 0
      return -1
    if x[4] <= x[2]: # y3 <= y2
      return -1
    if x[6] <= x[4]: # y4 <= y3
      return -1
    if x[12] <= x[6]: # hei <= y4
      return -1
    return 1
    
h0 = [0.0, 0.02, 0.01, -0.04, 0.02, 0.01, 0.03, 0.0, 0.0013, 0.0052, 0.0013, 0.04]
constraints = [g0, g1, g2]
h_opt = fmin_cobyla(objective, h0, constraints, rhoend=1e-6, maxfun=100, catol=1e-6)
results_write.write(str(h_opt))
results_write.close()
print h_opt, objective(h_opt), g0(h_opt)
