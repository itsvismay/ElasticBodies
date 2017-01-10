# inputs
# --name             name of the file to be outputed
# --minRadius        the minimum value of the radius
# --maxRadius        the maximum value of the radius
# --minHeight        the minimum value of the height
# --maxHeight        the maximum value of the height
# --minWidth         the minimum value of the width
# --maxWidth         the maximum value of the width
# --resolution       the total number of samples
# --loops            the fractional number of loops
# -a read raw data

# raw data is of this form:
# number of control points for radius
# number of control points for width
# number of control points for height
# min radius
# max radius
# min width
# max width
# min height
# max height
# resolution
# loops
# for each control point for radius
#   the control points y value
# for each control point for width
#   the control points y value
# for each control point for height
#   the control points y value

# there are 3 1D Bezier Curves to represnt this
# For radius
#   first control point is locked at 0.0
#   last control point is locked at 1.0
# For width
#   first control point is at 0.5
#   second control point is at 0.5
# For height
#   first control point is at 0.5
#   second control point is at 0.5

import sys, getopt, os, math
import numpy as np

name = "torsion_test.scad"
numberRadiusCtrlPts = 5
numberWidthCtrlPts = 5
numberHeightCtrlPts = 5
radiusCtrlPts = [0.0, 0.25, 0.5, 0.75, 1.0]
widthCtrlPts = [1.0, 1.0, 1.0, 1.0, 1.0]
heightCtrlPts = [0.2, 0.7, 0.95, 0.7, 0.2]

maxRadius = 30.0
minRadius = 0.0
minWidth = 1.0
maxWidth = 5.0
minHeight = 3.0
maxHeight = 10.0
resolution = 100
loops = 3.0

try:
  opts, args = getopt.getopt(sys.argv[1:], 'a', ["maxRadius=", "minRadius=", "minHeight=", "maxHeight=", "minWidth=", "maxWidth=", "resolution=", "loops="])
except getopt.GetoptError:
  print 'Error Bad Input'
  sys.exit(-2)
for opt, arg in opts:
  if opt == "--name":
    name = arg
  elif opt == "--maxRadius":
    maxRadius = float(arg)
  elif opt == "--minRadius":
    minRadius = float(arg)
  elif opt == "--minHeight":
    minHeight = float(arg)
  elif opt == "--maxHeight":
    maxHeight = float(arg)
  elif opt == "--minWidth":
    minWidth = float(arg)
  elif opt == "--maxWidth":
    maxWidth = float(arg)
  elif opt == "--resolution":
    resolution = int(float(arg))
  elif opt == "--loops":
    loops = float(arg)
  elif opt == "-a":
    # do input stuffs
    radiusCtrlPts = []
    widthCtrlPts = []
    heightCtrlPts = []

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
  print "NUM CTRL ::", len(ctrlPts)
  for i in range(0, resolution):
    t = float(i) / float(resolution)
    print "t :", t
    nt = 1.0 - t
    print 'nt :', nt
    n = numberCtrls - 1
    currentVal = 0.0
    for j in range(0, len(ctrlPts)):
      tmp = math.pow(t, j) * math.pow(nt, n-j)
      tmp = tmp * combination(n, j)
      currentVal = currentVal + tmp * ctrlPts[j]
    evaluatedValues.append(currentVal)
  return evaluatedValues

def calculateRadius(v):
  return (maxRadius - minRadius) * v + minRadius

def calculateHeight(v):
  return (maxHeight - minHeight) * v + minHeight

def calculateWidth(v):
  return (maxWidth - minWidth) * v + minWidth

evaluatedRadius = evaluateBez(radiusCtrlPts)
evaluatedHeight = evaluateBez(heightCtrlPts)
evaluatedWidth  = evaluateBez(widthCtrlPts)

Data = []

for i in range(0, len(evaluatedRadius)):
  theta = (2 * math.pi * loops / len(evaluatedRadius)) * i
  pointX = calculateRadius(evaluatedRadius[i]) * math.cos(theta)
  pointY = calculateRadius(evaluatedRadius[i]) * math.sin(theta)
  width = calculateWidth(evaluatedWidth[i])
  height = calculateHeight(evaluatedHeight[i])
  Data.append([pointX, pointY, width, height])

print "5 choose 2 ::", combination(5, 2)
print "9 choose 4 ::", combination(9, 4)
print "3 choose 1 ::", combination(3, 1)

print radiusCtrlPts
#print evaluatedRadius
#print evaluatedWidth
print Data

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
  file_write.write("[%7.8f, %8.8f, %7.8f, %7.8f],\n" % (row[0], row[1], row[2], row[3]))
file_write.write("];\n")
file_write.write("\n")
file_write.write("translate([0.0, 0.0, %6.8f]) {\n" % z)
file_write.write("\tlinear_extrude(height = %6.8f, center = false, convexity = 10, twist = 0) {\n" % (layerThk * 1.01))
file_write.write("\t\tunion() {\n")
file_write.write("\t\t\tfor( i=[0:len(points)-2] ){ \n")
file_write.write("\t\t\t\tx1 = points[i][0];\n")
file_write.write("\t\t\t\ty1 = points[i][1];\n")
file_write.write("\t\t\t\tx2 = points[i+1][0];\n")
file_write.write("\t\t\t\ty2 = points[i+1][1];\n")
file_write.write("\t\t\t\twidth = points[i+1][2];\n")
file_write.write("\t\t\t\ttranslateAndRotate(x1, y1, x2, y2, width);\n")
file_write.write("\t\t\t}\n")
file_write.write("\t\t}\n")
file_write.write("\t}\n")
file_write.write("}\n")
file_write.write("\n")
for row in Data:
  file_write.write("translate([%7.8f, %8.8f, 0.0]) circle(.1);\n" % (row[0], row[1]))
file_write.close()
