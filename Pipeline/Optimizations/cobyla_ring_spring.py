from scipy.optimize import fmin_cobyla
import sys, os, subprocess, numpy

# inputs / simulation constants
# --minThkInPlane            the minimum thickness in plane
# --maxThkInPlane            the maximum thickness in plane
# --minWidthInPlane          the minimum width in plane
# --maxWidthInPlane          the maximum width in plane
# --minThkOutOfPlane         the minimum thickness out of plane
# --maxThkOutOfPlane         the maximum thickness out of plane
# --minBeamWidth             the minimum width of the center beams
# --maxBeamWidth             the maximum width of the center beams
# --minHeight                the minimum height of the spring
# --maxHeignt                the maximum height of the spring
# --minSections              the minimum number of sections
# --maxSections              the maximum number of sections
# --minOverlap               the minimum amount of overlap
# --maxOverlap               the maximum amoung of overlap
# --minScaleX1               the minimum scaleX1 value
# --maxScaleX1               the maximum scaleX1 value
# --minScaleX2               the minimum scaleX2 value
# --maxScaleX2               the maximum scaleX2 value
# --minScaleY1               the minimum scaleY1 value
# --maxScaleY1               the maximum scaleY1 value
# --minScaleY2               the minimum scaleY2 value
# --maxScaleY2               the maximum scaleY2 value
# -a                         read everything from cmd arguments

minThkInPlane = 0.0
maxThkInPlane = 0.0
minWidthInPlane = 0.0
maxWidthInPlane = 0.0
minThkOutOfPlane = 0.0
maxThkOutOfPlane = 0.0
minBeamWidth = 0.0
maxBeamWidth = 0.0
minHeight = 0.0
maxHeight = 0.0
minSections = 0.0
maxSections = 0.0
minOverlap = 0.0
maxOverlap = 0.0
minScaleX1 = 0.0
maxScaleX1 = 0.0
minScaleX2 = 0.0
maxScaleX2 = 0.0
minScaleY1 = 0.0
maxScaleY1 = 0.0
minScaleY2 = 0.0
maxScaleY2 = 0.0

try:
    opts, args = getopt.getopt(sys.argv[1:], 'a', ["minThkInPlane=", "maxThkInPlane", "minWidthInPlane=", "maxWidthInPlane=", "minThkOutOfPlane=", "maxThkOutOfPlane=", "minBeamWidth=", "maxBeamWidth=", "minHeight=", "maxHeight=", "minSections=", "maxSections=", "minOverlap=", "maxOverlap=", "minScaleX1=", "maxScaleX1=", "minScaleX2=", "maxScaleX2=", "minScaleY1=", "maxScaleY1=", "minScaleY2=", "maxScaleY2="])
except getopt.GetoptError:
    print 'Error Bad Input'
    sys.exit(-2)
for opt, arg in opts:
    if opt == "--minThkInPlane":
		minThkInPlane = float(arg)
	elif opt == "--maxThkInPlane":
		maxThkInPlane = float(arg)
	elif opt == "--minWidthInPlane":
		minWidthInPlane = float(arg)
	elif opt == "--maxWidthInPlane":
		maxWidthInPlane = float(arg)
	elif opt == "--minThkOutOfPlane":
		minThkOutOfPlane = float(arg)
	elif opt == "--maxThkOutOfPlane":
		maxThkOutOfPlane = float(arg)
	elif opt == "--minBeamWidth":
		minBeamWidth = float(arg)
	elif opt == "--maxBeamWidth":
		maxBeamWidth = float(arg)
	elif opt == "--minHeight":
		minHeight = float(arg)
	elif opt == "--maxHeight":
		maxHeight = float(arg)
	elif opt == "--minSections":
		minSections = float(arg)
	elif opt == "--maxSections":
		maxSections = float(arg)
	elif opt == "--minOverlap":
		minOverlap = float(arg)
	elif opt == "--maxOverlap":
		maxOverlap = float(arg)
	elif opt == "--minScaleX1":
		minScaleX1 = float(arg)
	elif opt == "--maxScaleX1":
		maxScaleX1 = float(arg)
	elif opt == "--minScaleX2":
		minScaleX2 = float(arg)
	elif opt == "--maxScaleX2":
		maxScaleX2 = float(arg)
	elif opt == "--minScaleY1":
		minScaleY1 = float(arg)
	elif opt == "--maxScaleY1":
		maxScaleY1 = float(arg)
	elif opt == "--minScaleY2":
		minScaleY2 = float(arg)
	elif opt == "--maxScaleY2":
		maxScaleY2 = float(arg)
	elif opt == "-a":
        minThkInPlane = float(sys.argv[2])
        maxThkInPlane = float(sys.argv[3])
        minWidthInPlane = float(sys.argv[4])
        maxWidthInPlane = float(sys.argv[5])
        minThkOutOfPlane = float(sys.argv[6])
        maxThkOutOfPlane = float(sys.argv[7])
        minBeamWidth = float(sys.argv[8])
        maxBeamWidth = float(sys.argv[9])
        minHeight = float(sys.argv[10])
        maxHeight = float(sys.argv[11])
        minSections = int(float(sys.argv[12]))
        maxSections = int(float(sys.argv[13]))
        minOverlap = float(sys.argv[14])
        maxOverlap = float(sys.argv[15])
        minScaleX1 = float(sys.argv[16])
        maxScaleX1 = float(sys.argv[17])
        minScaleX2 = float(sys.argv[18])
        maxScaleX2 = float(sys.argv[19])
        minScaleY1 = float(sys.argv[20])
        maxScaleY1 = float(sys.argv[21])
        minScaleY2 = float(sys.argv[22])
        maxScaleY2 = float(sys.argv[23])
    else:
        print 'ERROR ::', arg, 'NOT RECOGNIZED'
        exit(10)

def objective(x):
    # to be implemented
    return 1.0

def g0(x):
    for i in range(0, 11):
        if x[i] < 0.0:
            return x[i]
    return 1.0

def g1(x):
    for i in range(0, 11):
        if x[i] > 1.0:
            return -x[i]
    return 1.0

h0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
constraints = [g0, g1]
h_opt = fmin_cobyla(objective, h0, constraints, rhoend=1e-6, maxfun=100, catol=1e-6)
print h_opt
