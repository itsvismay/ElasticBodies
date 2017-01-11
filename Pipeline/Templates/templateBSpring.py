# inputs
# --name             name of the file to be outputed
# --minRadius        the minimum value of the radius
# --maxRadius        the maximum value of the radius
# --height           the height of the spring
# --depth            the depth of the spring
# --minWidth         the minimum value of the width
# --maxWidth         the maximum value of the width
# --resolution       the total number of samples
# --loops            the fractional number of loops
# -a read raw data

# raw data is of this form:
# number of control points for radius
# number of control points for width
# min x
# max x
# min width
# max width
# height
# depth
# resolution
# for each control point for x
#   the control points y value
# for each control point for width
#   the control points y value

# there are 3 1D Bezier Curves to represnt this
# For x
#   first control point is locked at 0.5
#   last control point is locked at 0.5
# For width
#   first control point is locked at 0.0
#   last control point is locked at 0.0

import sys, getopt, os, math
import numpy as np

name = "bezier_test.scad"
numberXCtrlPts = 7
numberWidthCtrlPts = 7
xCtrlPts = [0.0, 0.25, 0.75, 0.2, 0.5, 0.9, 1.0]
widthCtrlPts = [0.6, 1.0, 0.7, 0.4, 0.2, 0.4, 0.5]

maxX = 40.0
minX = -40.0
minWidth = 1.0
maxWidth = 10.0
height = 40.0
depth = 10.0
resolution = 1000

try:
  opts, args = getopt.getopt(sys.argv[1:], 'a', ["maxX=", "minX=", "minWidth=", "maxWidth=", "height=", "resolution="])
except getopt.GetoptError:
  print 'Error Bad Input'
  sys.exit(-2)
for opt, arg in opts:
  if opt == "--name":
    name = arg
  elif opt == "--maxX":
    maxX = float(arg)
  elif opt == "--minX":
    minX = float(arg)
  elif opt == "--height":
    height = float(arg)
  elif opt == "--depth":
    depth = float(arg)
  elif opt == "--minWidth":
    minWidth = float(arg)
  elif opt == "--maxWidth":
    maxWidth = float(arg)
  elif opt == "--resolution":
    resolution = int(float(arg))
  elif opt == "-a":
    # do input stuffs
    numberXCtrlPts = int(float(sys.argv[3]))
    numberWidthCtrlPts = int(float(sys.argv[4]))
    minX = float(sys.argv[5])
    maxX = float(sys.argv[6])
    minWidth = float(sys.argv[7])
    maxWidth = float(sys.argv[8])
    height = float(sys.argv[9])
    resolution = int(float(sys.argv[10]))
    val = 10
    for i in range(0, len(numberXCtrlPts)):
      xCtrlPts.append(float(sys.argv[val+1]))
      val = val + 1
    for i in range(0, len(numberWidthCtrlPts)):
      widthCtrlPts.append(float(sys.argv[val+1]))
      val = val + 1

xCtrlPts.insert(0, 0.5)
xCtrlPts.append(0.5)
widthCtrlPts.insert(0, 0.0)
widthCtrlPts.append(0.0)
numberXCtrlPts += 2
numberWidthCtrlPts += 2

def factorial(v):
  val = 1;
  for i in range(1,v+1):
    val = val * i
  return val

def combination(n, r):
  return factorial(n) / (factorial(r) * factorial(n-r))

def evaluateBez(ctrlPts):
  evaluatedValues = []
  numberCtrls = len(ctrlPts)
  for i in range(0, resolution+1):
    t = float(i) / float(resolution)
    nt = 1.0 - t
    n = numberCtrls - 1
    currentVal = 0.0
    for j in range(0, len(ctrlPts)):
      tmp = math.pow(t, j) * math.pow(nt, n-j)
      tmp = tmp * combination(n, j)
      currentVal = currentVal + tmp * ctrlPts[j]
    evaluatedValues.append(currentVal)
  return evaluatedValues

def calculateX(v):
  return (maxX - minX) * v + minX

def calculateHeight(v):
  return height * v

def calculateWidth(v):
  return (maxWidth - minWidth) * v + minWidth

evaluatedX = evaluateBez(xCtrlPts)
evaluatedWidth  = evaluateBez(widthCtrlPts)

Data = []

for i in range(0, len(evaluatedX)):
  val = float(i) / float(len(evaluatedX))
  pointX = calculateX(evaluatedX[i])
  pointY = height * val
  width = calculateWidth(evaluatedWidth[i])
  Data.append([pointX, pointY, width])

layerThk = 0.4
z = 0

file_write = open(name, "w")
file_write.write("module drawBasicShape(x1, y1, x2, y2, width)\n")
file_write.write("{\n")
file_write.write("\tlength = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));\n")
file_write.write("\tif (width!=0.0) {\n")
file_write.write("\t\tsquare(size=[width, length], center = false);\n")
file_write.write("\t\ttranslate([width/2.0, 0, 0]) circle(d=width);\n")
file_write.write("\t\ttranslate([width/2.0, length, 0]) circle(d=width);\n")
file_write.write("\t}\n")
file_write.write("}\n")
file_write.write("\n")
file_write.write("module translateAndRotate(x1, y1, x2, y2, width)\n")
file_write.write("{\n")
file_write.write("\tif (width != 0.0) {\n")
file_write.write("\t\tangle = atan((y2-y1)/(x2-x1)) +270;\n")
file_write.write("\t\tif (x2>=x1) {\n")
file_write.write("\t\t\ttranslate([x1, y1, 0.0]) rotate([0.0, 0.0, angle]) translate([-width/2.0, 0.0, 0.0]) drawBasicShape(x1, y1, x2, y2, width);\n")
file_write.write("\t\t}\n")
file_write.write("\t\telse {\n")
file_write.write("\t\t\ttranslate([x1,y1, 0.0]) rotate([0.0, 0.0, angle+180]) translate([-width/2.0, 0.0, 0.0]) drawBasicShape(x1, y1, x2, y2, width);\n")
file_write.write("\t\t}\n")
file_write.write("\t}\n")
file_write.write("}\n")
file_write.write("\n")
file_write.write("$fn=50;\n")
file_write.write("\n")
file_write.write('points = [\n')
for row in Data:
  file_write.write("[%7.8f, %8.8f, %7.8f],\n" % (row[0], row[1], row[2]))
file_write.write("];\n")
file_write.write("\n")
file_write.write("translate([0.0, 0.0, %6.8f]) {\n" % z)
file_write.write("\tlinear_extrude(height = %6.8f, center = false, convexity = 10, twist = 0) {\n" % depth)
file_write.write("\t\tunion() {\n")
file_write.write("\t\t\tfor( i=[0:len(points)-2] ){ \n")
file_write.write("\t\t\t\tx1 = points[i][0];\n")
file_write.write("\t\t\t\ty1 = points[i][1];\n")
file_write.write("\t\t\t\tx2 = points[i+1][0];\n")
file_write.write("\t\t\t\ty2 = points[i+1][1];\n")
file_write.write("\t\t\t\twidth = points[i+1][2];\n")
file_write.write("\t\t\t\ttranslateAndRotate(x1, y1, x2, y2, width);\n")
file_write.write("\t\t\t}\n")
file_write.write("\t\t\ttranslate([%6.8f,%6.8f]) square([%6.8f, %6.8f]);\n" % (-(maxX / 2.0), -(minWidth), maxX, minWidth*2.0))
file_write.write("\t\t\ttranslate([%6.8f,%6.8f]) square([%6.8f, %6.8f]);\n" % (-(maxX / 2.0), height-(minWidth), maxX, minWidth*2.0))
file_write.write("\t\t}\n")
file_write.write("\t}\n")
file_write.write("}\n")
file_write.write("\n")
file_write.close()

print name
