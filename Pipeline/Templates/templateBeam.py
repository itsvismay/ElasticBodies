# inputs
# --name      name of file to be outputed
# --length    value of the length
# --width     value of the width
# --height    value of the height
# -a          read raw data

import sys, getopt, os

name = 'beam_test.scad'
length = 500
width = 100
height = 100

try:
  opts, args = getopt.getopt(sys.argv[1:], 'a', ["name=", "length=", "width=", "height="])
except getopt.GetoptError:
  print 'Error Bad Input'
  sys.exit(-2)
for opt, arg in opts:
  if opt == "--name":
    name = arg
  elif opt == "--length":
    length = float(arg)
  elif opt == "--width":
    width = float(arg)
  elif opt == "--height":
    height = float(arg)
  elif opt == "-a":
    name = sys.argv[2]
    length = float(sys.argv[3])
    width = float(sys.argv[4])
    height = float(sys.argv[5])

file_write = open(name, 'w')
file_write.write("$fn = 200;\n")
file_write.write("\n")
file_write.write("length = "+str(length)+";\n")
file_write.write("width = "+str(width)+";\n")
file_write.write("height = "+str(height)+";\n")
file_write.write("\n")
file_write.write("translate([0, 0, -height]) cube([length, width, height], circle = false);\n")
file_write.close()

print name
