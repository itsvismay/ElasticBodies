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
  print 'Error Bad Input'
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
    # my hypothesis
    deltaE = self.startE - self.endE
    if deltaE < -1.0:
      print 'ERROR ' + str(deltaE) + '\n'
    deltaD = np.sqrt((self.endX - self.startX)**2.0 + (self.endY - self.startY)**2)
    if deltaD == 0.0:
      return 0.0
    filArea = np.pi * (0.5 * filimentDiameter)**2.0
    heightFactor = layerHeight**2.0 * (np.pi / 4.0 - 1.0)
    return ( (deltaE / deltaD) * filArea - heightFactor ) / layerHeight
    #return 0.5

  def collidesWithY(self, y):
    fourCorners = self.getFourCorners() # can improve performance by caching
    maxY = -1000
    minY = 1000
    for i in range(len(fourCorners)):
      if fourCorners[i][1] > maxY:
        maxY = fourCorners[i][1]
      if fourCorners[i][1] < minY:
        minY = fourCorners[i][1]
    return y >= minY and y <= maxY

  def collidesWithPoint(self, y, x):
    #fourCorners = self.getFourCorners() # can improve performance by caching
    distVec = np.array([self.endX - self.startX, self.endY - self.startY])
    pointVec = np.array([x - self.startX, y - self.startY])
    moveVec = np.cross(distVec, np.array([0.0, 0.0, 1.0]))
    # projDistVec = (pointVec dot norm(distVec)
    # projMovVec = ((projDistVec - pointVec) dot norm(moveVec))
    # to be implemented
    return False

  def getFourCorners(self):
    fourCorners = []
    halfStartWidth = self.startW / 2.0
    halfEndWidth = self.endW / 2.0
    distVec = np.array([self.endX - self.startX, self.endY - self.startY, 0.0])
    moveVec = np.cross(distVec, np.array([0.0, 0.0, 1.0]))
    fourCorners.append([self.startX + moveVec[0], self.startY + moveVec[1]])
    fourCorners.append([self.startX - moveVec[0], self.startY - moveVec[1]])
    fourCorners.append([self.endX - moveVec[0], self.endY - moveVec[1]])
    fourCorners.append([self.endX + moveVec[0], self.endY + moveVec[1]])
    return fourCorners

class Row:
  extrusionsInRow = []
  rowData = []
  rowY = 0.0
  dim = 0.0
  numRow = 0
  layer = -1

  def __init__(self, extrusions, y, d, numR, startX, lay):
    self.rowY = y
    self.dim = d
    self.numRow = numR
    self.layer = lay
    getExtrusionsInRow()

  def getExtrusionsInRow(self, extrusions, startX):
    for ex in extrusions:
      if (ex.layer == self.layer):
        if (ex.collidesWithY(self.rowY) == True):
          rowData.append(ex)

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
          #print str(tokenE), ' :: ', str(prevE)
          newExtrusion = Extrusion(tokenX, tokenY, tokenE, prevF, prevX, prevY, prevE, prevF, layer)
          #print str(tokenE), ' :: ', str(prevE)
          extrusions.append([layer, newExtrusion])
          newData.append([layer, tokenX, tokenY, newExtrusion.startW]) # replace with the extrusion width soon
        elif (tokenX > -1.0 and tokenY > -1.0):
          newData.append([layer, tokenX, tokenY, 0.0])
      #if (tokenF > -1.0):
      #  prevF = tokenF
      if (tokenE > -1.0):
        prevE = tokenE
      if (tokenX > -1.0):
        prevX = tokenX
      if (tokenY > -1.0):
        prevY = tokenY
      if (tokenZ > -1.0):
        prevX = -1.0
        prevY = -1.0
        #prevE = -1.0
        layer += 1

  return newData, extrusions, count, layer

def getVoxelDimensions(extrusions):
  minX = 10000
  minY = 10000
  maxX = -10000
  maxY = -10000
  minW = 10000

  for ex in extrusions:
    if (minX > ex.startX):
      minX = ex.startX
    if (minX > ex.endX):
      minX = ex.endX
    if (minY > ex.startY):
      minY = ex.startY
    if (minY > ex.endY):
      minY = ex.endY
    if (maxX < ex.startX):
      maxX = ex.startX
    if (maxX < ex.endX):
      maxX = ex.endX
    if (maxY < ex.startY):
      maxY = ex.startY
    if (maxY < ex.endY):
      maxY = ex.endY
    if (minW > ex.startW):
      minW = ex.startW

  dim = minW / 3.0
  return dim, minX, maxX, minY, maxY

print vName
Data_old, numLines_old = loadData(vName); print numLines_old
Data, extrusions, numLines, numLayers = cleanData(Data_old, numLines_old)
print numLines, numLayers

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
  filename = "CCgcode_layer" + str(layer) + ".scad"
  file_write = open(filename, "w")
  file_write.write('points' + str(layer) + ' = [\n')
  for row in layerData:
    file_write.write("[%7.3f, %8.3f, %2.6f],\n" % (row[1], row[2], row[3]))
  file_write.write("];\n")
  z += layerThk
  file_write.write("translate([0.0, 0.0, %6.3f]) {\n" % (z))
  file_write.write("\tlinear_extrude(height = %6.3f, center = false, convexity = 10, twist = 0) {\n" % (layerThk))
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

filename = "CCgcode_for_loop.scad"
file_write = open(filename, "w")
for layer in runLayers:
  file_write.write("include <CCgcode_layer%d.scad>;\n" % (layer))
file_write.close()
