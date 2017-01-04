# inputs
# --name + file name

import sys, getopt, os
import numpy as np

vName = 'TestImage.gcode'
filimentDiameter = 1.75
layerHeight = 0.4

try:
  opts, args = getopt.getopt(sys.argv[1:], 'f', ["name="])
except getopt.GetoptError:
  #print 'Error Bad Input'
  sys.exit(-2)
for opt, arg in opts:
  if opt == "--name":
    vName = arg

class Extrusion:
  startX = 0.0 # xpos
  startY = 0.0 # ypos
  startE = 0.0 # extrusion length
  startF = 0.0 # feed rate
  startW = 0.0 # width
  endX = 0.0
  endY = 0.0
  endE = 0.0
  endF = 0.0
  endW = 0.0
  layer = -1

  def __init__(self, sX, sY, sE, sF, eX, eY, eE, eF, lay):
    self.startX = sX
    self.startY = sY
    self.startE = sE
    self.startF = sF
    self.endX = eX
    self.endY = eY
    self.endE = eE
    self.endF = eF
    #print str(self.startE), ' :2: ', str(self.endE)
    self.layer = lay
    self.startW = self.calculateWidth(sX, sY, sF, sE, eE)
    self.endW = self.calculateWidth(eX, eY, eF, sE, eE)

  def calculateWidth(self, x, y, f, e, pE):
    #return 0.5
    deltaE = self.startE - self.endE
    #if deltaE < -1.0:
    #  print 'ERROR ' + str(deltaE) + '\n'
    deltaD = np.sqrt((self.endX - self.startX)**2.0 + (self.endY - self.startY)**2)
    if deltaD == 0.0:
      return 0.0
    filArea = np.pi * (0.5 * filimentDiameter)**2.0
    heightFactor = layerHeight**2.0 * (np.pi / 4.0 - 1.0)
    return ( (deltaE / deltaD) * filArea - heightFactor ) / layerHeight

def loadData(fileName):
  fr = open(fileName)
  Data = []; count = 0
  for line in fr.readlines():
    count += 1
    lineArr = []
    curLine = line.strip().split(' ')
    for i in range(len(curLine)):
      lineArr.append(curLine[i])
    Data.append(lineArr)
  return Data, count

def cleanData(rawData, numLines):

  newData = []
  extrusions = []
  count = 0;
  layer = -1;
  prevX = -1.0
  prevY = -1.0
  prevE = -1.0
  prevF = -1.0

  for i in range(numLines):
    if ("G1" in rawData[i][0]):
      if (";" in rawData[i]):
        continue
      tokenX = -1.0
      tokenY = -1.0
      tokenE = -1.0
      tokenF = -1.0
      tokenZ = -1.0
      for j in range(len(rawData[i])):
        if("Z" in rawData[i][j]):
          tokenZ = float(rawData[i][j].strip('Z'))
        if("X" in rawData[i][j]):
          tokenX = float(rawData[i][j].strip('X'))
        if("Y" in rawData[i][j]):
          tokenY = float(rawData[i][j].strip("Y"))
        if("F" in rawData[i][j]):
          tokenF = float(rawData[i][j].strip("F"))
        if("E" in rawData[i][j]):
          tokenE = float(rawData[i][j].strip("E"))
      if (tokenF > -1.0):
        prevF = tokenF
      if (not (tokenX == prevX and tokenY == prevY and tokenE == prevE) ):
        if (tokenX > -1.0 and tokenY > -1.0 and tokenE > -1.0 and tokenZ == -1.0):
          # create extrusion
          newExtrusion = Extrusion(tokenX, tokenY, tokenE, prevF, prevX, prevY, prevE, prevF, layer)
          newData.append([layer, tokenX, tokenY, newExtrusion.startW]) # replace with the extrusion width soon
        elif (tokenX > -1.0 and tokenY > -1.0):
          newData.append([layer, tokenX, tokenY, 0.0])
      if (tokenE > -1.0):
        prevE = tokenE
      if (tokenX > -1.0):
        prevX = tokenX
      if (tokenY > -1.0):
        prevY = tokenY
      if (tokenZ > -1.0):
        prevX = -1.0
        prevY = -1.0
        layer += 1

  return newData, extrusions, count, layer

#print vName
Data_old, numLines_old = loadData(vName); #print numLines_old
Data, extrusions, numLines, numLayers = cleanData(Data_old, numLines_old)
#print numLines, numLayers

#-------------------------------------------------
# Code for writing the layer files into scad
#-------------------------------------------------
z = 0
layerThk = 0.4
indices = np.array([row[0] for row in Data])
DataNP = np.array(Data)
runLayers = range(numLayers+1)
for layer in runLayers:
  layerData = DataNP[indices == layer]
  filename = vName[:-6]+"_layer_" + str(layer) + ".scad"
  file_write = open(filename, "w")
  # write needed modules
  file_write.write("module drawBasicShape(x1, y1, x2, y2, width)\n")
  file_write.write("{\n")
  file_write.write("\tlength = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));\n")
  file_write.write("\tif (width!=0.0) {\n")
  file_write.write("\t\tsquare(size = [width, length], center = false);\n")
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
  if layer == 0:
    file_write.write("\t\t\ttranslate([x1, y1, 0.0]) rotate([0.0, 0.0, angle]) translate([-width/2.0, 0.0, 0.0]) drawBasicShape(x1, y1, x2, y2, width);\n")
  else:
    file_write.write("\t\t\ttranslate([x1, y1, 0.0]) rotate([0.0, 0.0, angle]) translate([-width/2.0, 0.0, 0.0]) scale([1.0, 1.0, 1.0]) drawBasicShape(x1, y1, x2, y2, width);\n")
  file_write.write("\t\t}\n")
  file_write.write("\t\telse {\n")
  if layer == 0:
    file_write.write("\t\t\ttranslate([x1, y1, 0.0]) rotate([0.0, 0.0, angle+180]) translate([-width/2.0, 0.0, 0.0]) drawBasicShape(x1, y1, x2, y2, width);")
  else:
    file_write.write("\t\t\ttranslate([x1, y1, 0.0]) rotate([0.0, 0.0, angle+180]) translate([-width/2.0, 0.0, 0.0]) scale([1.0, 1.0, 1.0]) drawBasicShape(x1, y1, x2, y2, width);")
  file_write.write("\t\t}\n")
  file_write.write("\t}\n")
  file_write.write("}\n")
  file_write.write("\n")
  file_write.write("$fn=50;\n")
  file_write.write("\n")

  # write actual data
  file_write.write('points' + str(layer) + ' = [\n')
  for row in layerData:
    file_write.write("[%7.8f, %8.8f, %2.8f],\n" % (row[1], row[2], row[3]))
  file_write.write("];\n")
  file_write.write("\n")
  z += layerThk
  #file_write.write("translate([0.0, 0.0, %6.3f]) {\n" % (z - (layerThk * .005)))
  file_write.write("translate([0.0, 0.0, %6.8f]) {\n" % z)
  if layer == 0:
    file_write.write("\tlinear_extrude(height = %6.8f, center = false, convexity = 10, twist = 0) {\n" % (layerThk))
  else:
    file_write.write("\tlinear_extrude(height = %6.8f, center = false, convexity = 10, twist = 0) {\n" % (layerThk * 1.01))
  file_write.write("\t\tunion() {\n")
  file_write.write("\t\t\tfor( i=[0:len(points" + str(layer) + ")-2] ){ \n")
  file_write.write("\t\t\t\tx1 = points" + str(layer) + "[i][0];\n")
  file_write.write("\t\t\t\ty1 = points" + str(layer) + "[i][1];\n")
  file_write.write("\t\t\t\tx2 = points" + str(layer) + "[i+1][0];\n")
  file_write.write("\t\t\t\ty2 = points" + str(layer) + "[i+1][1];\n")
  file_write.write("\t\t\t\twidth = points" + str(layer) + "[i+1][2];\n")
  file_write.write("\t\t\t\ttranslateAndRotate(x1, y1, x2, y2, width);\n")
  file_write.write("\t\t\t}\n")
  file_write.write("\t\t}\n")
  file_write.write("\t}\n")
  file_write.write("}\n")
  file_write.close()

#filename = "CCgcode_for_loop.scad"
#file_write = open(filename, "w")
#for layer in runLayers:
#  file_write.write("include <"+vName[:-6]+"_layer_%d.scad>;\n" % (layer))
#file_write.close()

print numLayers

#sys.exit(numLayers)
