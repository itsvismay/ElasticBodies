from scipy.optimize import fmin_cobyla
import sys, os, subprocess, numpy

# inputs / simulation constants
# --minInThk      minimum in plane thickness
# --maxInThk      maximum in plane thickness
# --minInHei      minimum in height
# --maxInHei      maximum in height
# --minSections   minimum number of sections
# --maxSections   maximum number of sections
# --minOutThk     minimum out of plane thickness
# --maxOutThk     maximum out of plane thickness
# --minWidthIn    minimum width in plane
# --maxWidthIn    maximum width in plane
# --minOverlap    minimum overlap
# --maxOverlap    maximum overlap
# -a                 read everything from cmd arguments

minInThk = 0.0
maxInThk = 0.0
minInHei = 0.0
maxInHei = 0.0
minSections = 0.0
maxSections = 0.0
minOutThk = 0.0
maxOutThk = 0.0
minWidthIn = 0.0
maxWidthIn = 0.0
minOverlap = 0.0
maxOverlap = 0.0

try:
    opts, args = getopt.getopt(sys.argv[1:], 'a', ["minInThk=", "maxInThk=", "minInHei=", "maxInHei=", "minSections=", "maxSections=", "minOutThk=", "maxOutThk", "minWidthIn=", "maxWidthIn=", "minOverlap=", "maxOverlap="])
except getopt.GetoptError:
    print 'Error Bad Input'
    sys.exit(-2)
for opt, arg in opts:
    if opt == "--minInThk":
        minInThk = float(arg)
    elif opt == "--maxInThk":
        maxInThk = float(arg)
    elif opt == "--minInHei":
        minInHei = float(arg)
    elif opt == "--maxInHei":
        maxInHei = float(arg)
    elif opt == "--minSections":
        minSections = float(arg)
    elif opt == "--maxSections":
        maxSections = float(arg)
    elif opt == "--minOutThk":
        minOutThk = float(arg)
    elif opt == "--maxOutThk":
        maxOutThk = float(arg)
    elif opt == "--minWidthIn":
        minWidthIn = float(arg)
    elif opt == "--maxWidthIn":
        maxWidthIn = float(arg)
    elif opt == "--minOverlap":
        minOverlap = float(arg)
    elif opt == "--maxOverlap":
        maxOverlap = float(arg)
    elif opt == "-a":
        minInThk = float(sys.argv[2])
        maxInThk = float(sys.argv[3])
        minInHei = float(sys.argv[4])
        maxInHei = float(sys.argv[5])
        minSections = float(sys.argv[6])
        maxSections = float(sys.argv[7])
        minOutThk = float(sys.argv[8])
        maxOutThk = float(sys.argv[9])
        minWidthIn = float(sys.argv[10])
        maxWidthIn = float(sys.argv[11])
        minOverlap = float(sys.argv[12])
        maxOverlap = float(sys.argv[13])
    else:
        print 'ERROR ::', arg, 'NOT RECOGNIZED'
        exit(10)

def objective(x):
    # to be implemented
    return 1.0

def g0(x):
    for i in range(0, 6):
        if x[i] < 0.0:
            return x[i]
    return 1.0

def g1(x):
    for i in range(0, 6):
        if x[i] > 1.0:
            return -x[i]
    return 1.0

h0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
constraints = [g0, g1]
h_opt = fmin_cobyla(objective, h0, constraints, rhoend=1e-6, maxfun=100, catol=1e-6)
print h_opt
