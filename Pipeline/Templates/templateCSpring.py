# inputs
# --name          name of the file to be outputed
### to be implemented later maybe ###
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
#####################################
# --inThk         thickness in plane
# --inHei         height in plane
# --sections      number of sections
# --outThk        thickness out of plane
# --widthIn       width in plane
# -a              read all vars as arguments

import sys, getopt, os, math
import numpy as np

name = 'curve_spring.scad'

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

inThk = 1.2
inHei = 40.0
sections = 5.0
outThk = 10.0
widthIn = 10.0

try:
  opts, args = getopt.getopt(sys.argv[1:], 'a', ["name=", "minInThk=", "maxInThk=", "minInHei=", "maxInHei=", "minSections=", "maxSections=", "minOutThk=", "maxOutThk=", "minWidthIn=", "maxWidthIn=", "inThk=", "inHei", "sections=", "outThk=", "widthIn="])
except getopt.GetoptError:
  print 'Error Bad Input'
  sys.exit(-2)
for opt, arg in opts:
  if opt == "--name":
    name = arg
  elif opt == "--minInThk":
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
  elif opt == "--inThk":
    inThk = float(arg)
  elif opt == "--inHei":
    inHei = float(arg)
  elif opt == "--sections":
    sections = float(arg)
  elif opt == "--outThk":
    outThk = float(arg)
  elif opt == "--widthIn":
    widthIn = float(arg)
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
    inThk = calculateInThk(float(sys.argv[12]))
    inHei = calculateInHei(float(sys.argv[13]))
    sections = calculateSections(float(sys.argv[14]))
    outThk = calculateOutThk(float(sys.argv[15]))
    widthIn = calculateWidthIn(float(sys.argv[16]))

def calculateInThk(v):
  return (maxInThk - minInThk) * v + minInThk

def calculateInHei(v):
  return (maxInHei - minInHei) * v + minInHei

def calculateSections(v):
  return (maxSections - minSections) * v + minSections

def calculateWidthIn(v):
  return (maxWidthIn - minWidthIn) * v + minWidthIn

def calculateOutThk(v):
  return (maxOutThk - minOutThk) * v + minOutThk

file_write = open(name, 'w')
file_write.write("thkInPlane = %3.8f;\n" % inThk)
file_write.write("widthInPlane = %3.8f;\n" % widthIn)
file_write.write("heightInPlane = %3.8f;\n" % inHei)
file_write.write("numberOfSectionsSet = %3.8f;\n"% sections)
file_write.write("thkOutOfPlane = %3.8f;\n" % outThk)
file_write.write("numberOfSections = floor(numberOfSectionsSet);\n")
file_write.write("radius = 0.5 * (heightInPlane/numberOfSections + thkInPlane/2);\n")
file_write.write("$fn=75;\n")
file_write.write("\n")
file_write.write("spring();\n")
file_write.write("\n")
file_write.write("module halfCylinder(rad, thk)\n")
file_write.write("{\n")
file_write.write("\tdifference() {\n")
file_write.write("\t\tdifference() {\n")
file_write.write("\t\t\tscale([1,1,1]) cylinder(r=radius, h=thk, center=true);\n")
file_write.write("\t\t\tscale([0.9,1,1]) cylinder(r=radius-thkInPlane, h=thk+2, center=true);\n")
file_write.write("\t\t};\n")
file_write.write("\ttranslate([-rad,0,0]) cube([2*rad,2*rad,thk+2], center=true);\n")
file_write.write("\t};\n")
file_write.write("}\n")
file_write.write("\n")
file_write.write("module spring() {\n")
file_write.write("cube([widthInPlane-radius, thkInPlane, thkOutOfPlane], center=false);\n")
file_write.write("for (i=[1:numberOfSections]) {\n")
file_write.write("\techo(i);\n")
file_write.write("\tif (i%2 != 0) {\n")
file_write.write("\t\techo(\"in the first if-loop\");\n")
file_write.write("\t\ttranslate([widthInPlane-radius,2*(i-1)*radius+radius-(i-1)*thkInPlane,thkOutOfPlane/2])\n")
file_write.write("\t\thalfCylinder(radius, thkOutOfPlane);\n")
file_write.write("\t\ttranslate([radius,2*i*radius-i*thkInPlane,0])\n")
file_write.write("\t\tcube([widthInPlane-2*radius,thkInPlane,thkOutOfPlane], center=false);\n")
file_write.write("\t}\n")
file_write.write("\telse {\n")
file_write.write("\t\techo(\"in the else-loop\");\n")
file_write.write("\t\ttranslate([radius,(2*i-1)*radius-(i-1)*thkInPlane,thkOutOfPlane/2])\n")
file_write.write("\t\trotate([0,0,180])\n")
file_write.write("\t\thalfCylinder(radius, thkOutOfPlane);\n")
file_write.write("\t\ttranslate([radius,2*i*radius-i*thkInPlane,0])\n")
file_write.write("\t\tcube([widthInPlane-2*radius,thkInPlane,thkOutOfPlane], center=false);\n")
file_write.write("\t}\n")
file_write.write("}\n")
file_write.write("if (numberOfSections%2 != 0) {\n")
file_write.write("translate([0,numberOfSections*(2*radius-thkInPlane),0])\n")
file_write.write("cube([radius, thkInPlane, thkOutOfPlane], center=false);\n")
file_write.write("}\n")
file_write.write("else {\n")
file_write.write("translate([widthInPlane-radius,numberOfSections*(2*radius-thkInPlane),0])\n")
file_write.write("\tcube([radius, thkInPlane, thkOutOfPlane], center=false);\n")
file_write.write("}\n")
file_write.write("}\n")
file_write.write("\n")
file_write.close()
