from scipy.optimize import fmin_cobyla
import sys, os, subprocess, numpy

# inputs / simulation constants
# --minRadius        the minimum value of the radius
# --maxRadius        the maximum value of the radius
# --height           the height of the spring
# --depth            the depth of the spring
# --minWidth         the minimum value of the width
# --maxWidth         the maximum value of the width
# --numberOfCtrlPts  the number of control points for radius and width
# -a                 read everything from cmd arguments

minRadius = 0.0
maxRadius = 0.0
height = 0.0
depth = 0.0
minWidth = 0.0
maxWidth = 0.0
resolution = 100
numberOfCtrlPts = 7

try:
    opts, args = getopt.getopt(sys.argv[1:], 'a', ["minRadius=", "maxRadius=", "height=", "depth=", "minWidth=", "maxWidth=", "resolution=", "numberOfCtrlPts="])
except getopt.GetoptError:
    print 'Error Bad Input'
    sys.exit(-2)
for opt, arg in opts:
    if opt == "--minRadius":
        minRadius = float(arg)
    elif opt == "--maxRadius":
        maxRadius = float(arg)
    elif opt == "--height":
        height = float(arg)
    elif opt == "--depth":
        depth = float(arg)
    elif opt == "--minWidth":
        minWidth = float(arg)
    elif opt == "--maxWidth":
        maxWidth = float(arg)
    elif opt == "--numberOfCtrlPts":
        numberOfCtrlPts = int(float(arg))
    elif opt == "--resolution":
        resolution = int(float(arg))
    elif opt == "-a":
        minRadius = float(sys.argv[2])
        maxRadius = float(sys.argv[3])
        height = float(sys.argv[4])
        depth = float(sys.argv[5])
        minWidth = float(sys.argv[6])
        maxWidth = float(sys.argv[7])
        resolution = int(float(sys.argv[8]))
        numberOfCtrlPts = int(float(sys.argv[9]))
    else:
        print 'ERROR ::', arg, 'NOT RECOGNIZED'
        exit(10)

def objective(x):
    # to be implemented
    return 1.0

def g0(x):
    for i in range(0, numberOfCtrlPts*2):
        if x[i] < 0.0:
            return x[i]
    return 1.0

def g1(x):
    for i in range(0, numberOfCtrlPts*2):
        if x[i] > 1.0:
            return -x[i]
    return 1.0

h0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
constraints = [g0, g1]
h_opt = fmin_cobyla(objective, h0, constraints, rhoend=1e-6, maxfun=100, catol=1e-6)
print h_opt
