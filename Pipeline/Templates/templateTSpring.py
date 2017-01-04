# inputs
# --name             name of the file to be outputed
# --ctrlPts          number of control points followed by the points x1 y1 x2 y2
# --width            width for each control point
# --height           height for each control point
# --displacement     displacement from center for each control point
# --resolution       the resolution (number of draws between control points)
# --dTheta           the change in theta per control point
# --complex          if the curve is rendererd by bspline or normal
# -a                 read raw data

import sys, getopt, os, math
import numpy as np

name = 'torsion_test.scad'
dataFile = 'torsionData.txt'
numCtrlPts = 20
ctrlPts = []
width = [1.0]
height = [1.0]
displacement = [0.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
complexity = 0 # 0 = simple 1 = bspline
resolution = 5
dTheta = 60

try:
  opts, args = getopt.getopt(sys.argv[1:], 'a', ["name=", "data=", "complex=", "resolution=", "displacement=", "dtheta="])
except getopt.GetoptError:
  print 'Error Bad Input'
  sys.exit(-2)
for opt, arg in opts:
  if opt == "--name":
    name = arg
  elif opt == "--data":
    dataFile = arg
  elif opt == "--resolution":
    resolution = float(arg)
  #elif opt == "--displacement":
  #  displacement = []
  elif opt == "--complex":
    complexity = int(arg)
  elif opt == "--dTheta":
    dTheta = float(arg)

def loadData(fileName):
  #fr = open(filename)
  count = 0;
  #for line in fr.readlines():
  #  curLine = line.strip().split(' ')
  #  if count == 0:
  #    # read in number of points
  #    temp = 0
  #  elif count == 1:
  #    # read in displacement
  #    temp = 0
  #  elif count == 2:
  #    # read in width
  #    temp = 0
  #  elif count == 3:
  #    # read in height
  #    temp = 0
  #  count += 1
  currentR = 0.0
  for i in range(numCtrlPts):
    currentD = displacement[i]
    currentR += currentD
    ctrlPts.append([currentR * math.cos(math.radians(dTheta * i)), currentR * math.sin(math.radians(dTheta * i))])
    print ctrlPts[i]

def createKnot(ctrlPts, degree, numCtrlPts):
  #val = len(ctrlPts) + degree + 1
  val = numCtrlPts + degree + 1
  #print 'Val', val
  siz = float(val - 1)
  knot = []
  for i in range(val):
    knot.append(float(i) / siz)
  return knot

def basis(i, j, knot, t):
  if j == 0 and t < knot[i+1] and t >= knot[i]:
    return 1
  if j == 0:
    return 0
  return ((t - knot[i]) / (knot[i+j] - knot[i])) * basis(i,j-1,knot,t) + ((knot[i+j+1] - t) / (knot[i+j+1] - knot[i+1])) * basis(i+1,j-1,knot,t)

def simple():
  Data = []
  #Data.append([0.0, 0.0, 1.0])
  #totalTheta = dTheta * numCtrlPts
  #ddTheta = dTheta / resolution
  #currentTheta = 0
  #for i in range(numCtrlPts):
  #  for j in range(resolution):
      
  # to be implemented

def bspline():
  Data = []
  #Data.append([0.0, 0.0, 1.0])
  knot = createKnot(ctrlPts, 3, numCtrlPts)
  print 'Knot'
  for i in range(len(knot)):
    print knot[i]
  samples = numCtrlPts * (resolution)
  print 'Samples', samples
  totalTheta = dTheta * numCtrlPts;
  print 'TotalTheta', totalTheta
  for i in range(samples):
    currentJ = 3
    currentI = math.floor(float(i) / float(resolution))
    currentT = float(i) / float(samples)
    #basisVal = basis(int(currentI), int(currentJ), knot, currentT)
    pointX = 0.0
    pointY = 0.0
    for j in range(numCtrlPts):
      basisVal = basis(int(j), int(currentJ), knot, currentT)
      pointX += ctrlPts[j][0] * basisVal
      pointY += ctrlPts[j][1] * basisVal
    Data.append([pointX, pointY, 1.0])
  return Data

# processing stuff goes here
loadData("data.txt")
data = bspline()
for i in range(len(data)):
  print data[i]
Data = np.array(data)

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
  file_write.write("[%7.8f, %8.8f, %2.8f],\n" % (row[0], row[1], row[2]))
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
for row in ctrlPts:
  file_write.write("translate([%7.8f, %8.8f, 0.0]) circle(.1);\n" % (row[0], row[1]))
file_write.close()
