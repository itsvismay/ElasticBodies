from scipy.optimize import fmin_cobyla
import sys, os, subprocess, numpy, getopt

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
# --maxWidthIn   maximum width in plane
# -a                 read everything from cmd arguments

minInThk = 0.4
maxInThk = 4.0
minInHei = 12.0
maxInHei = 38.0
minSections = 4.0
maxSections = 10.0
minOutThk = 9.0
maxOutThk = 12.0
minWidthIn = 9.0
maxWidthIn = 12.0

baseName = 'Individual'
counter = 0
TEST = 3
fileName = 'optimizeTest.txt'

try:
    opts, args = getopt.getopt(sys.argv[1:], 'a', ["minInThk=", "maxInThk=", "minInHei=", "maxInHei=", "minSections=", "maxSections=", "minOutThk=", "maxOutThk", "minWidthIn=", "maxWidthIn="])
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
    else:
        print 'ERROR ::', arg, 'NOT RECOGNIZED'
        exit(10)

def objective(x):
    file_write = open(fileName, 'w')
    file_write.write(fileName + ".scad "+str(minInThk)+" "+str(maxInThk)+" "+str(minInHei)+" "+str(maxInHei)+" "+str(minSections)+" "+str(maxSections)+" "+str(minOutThk)+" "+str(maxOutThk)+" "+str(minWidthIn)+" "+str(maxWidthIn)+" "+str(x[0])+" "+str(x[1])+" "+str(x[2])+" "+str(x[3])+" "+str(x[4]))
    file_write.close()

    subprocess.check_output(['python', 'pipeline.py', '--template', 'Templates/templateCSpring.py', '--create', 'optimizeTest.txt', '--sConfig', 'slic3rConfig.ini', '--preped', baseName+str(TEST), '-s', '-c'])

    
    return 1.0

def g0(x):
    for i in range(0, 5):
        if x[i] < 0.0:
            return x[i]
    return 1.0

def g1(x):
    for i in range(0, 5):
        if x[i] > 1.0:
            return -x[i]
    return 1.0

h0 = [0.5, 0.5, 0.5, 0.5, 0.5]
constraints = [g0, g1]
h_opt = fmin_cobyla(objective, h0, constraints, rhoend=1e-6, maxfun=100, catol=1e-6)
print h_opt
