# inputs
# --name    name of the file to be outputed
# --x3      value of the x coord of c3
# --y3      value of the y coord of c3
# -a        read raw data

import sys, getopt, os

name = 'beam_test.scad'
x3 = -40.0
y3 = 20.0

try:
  opts, args = getopt.getopt(sys.argv[1:], 'a', ["name=", "x3=", "y3="])
except getopt.GetoptError:
  print 'Error Bad Input'
  sys.exit(-2)
for opt, arg in opts:
  if opt == "--name":
    name = arg
  elif opt == "--x3":
    x3 = float(arg)
  elif opt == "--y3":
    y3 = float(arg)
  elif opt == "-a":
    name = sys.argv[2]
    x3 = float(sys.argv[3])
    y3 = float(sys.argv[4])

file_write = open(name, 'w')
file_write.write("use <BezierScad.scad>;\n")
file_write.write("bezierSpring(x3="+str(x3)+",y3="+str(y3)+");\n");
file_write.write("$fn=40;\n")
file_write.write("\n")
file_write.write("module bezierSpring(x3,y3) {\n")
file_write.write("\tx1=0.0;\n")
file_write.write("\tx2=20.0;\n")
file_write.write("\ty2=10.0;\n")
file_write.write("\tx4=10.0;\n")
file_write.write("\ty4=30.0;\n")
file_write.write("\tx5=0.0;\n")
file_write.write("\twidth2=1*1.30;\n")
file_write.write("\twidth3=4*1.30;\n")
file_write.write("\twidth4=1*1.30;\n")
file_write.write("\n")
file_write.write("\tthk=1.30;\n")
file_write.write("\tdia1=4.30;\n")
file_write.write("\tdia2=dia1-2*thk;\n")
file_write.write("\theight=40.0;\n")
file_write.write("\tunion() {\n")
file_write.write("\t\tbezierCurve();\n")
file_write.write("\t\tlowerPart();\n")
file_write.write("\t\ttranslate([0,height+dia1-thk/2.0]) upperPart();\n")
file_write.write("\t}\n")
file_write.write("\n")
file_write.write("\tmodule bezierCurve() {\n")
file_write.write("\t\tlinear_extrude(height=10)\n")
file_write.write("\t\tBezLine([\n")
file_write.write("\t\t\t[x1,-thk/2],\n")
file_write.write("\t\t\t[x2, y2],\n")
file_write.write("\t\t\t[x3,y3],\n")
file_write.write("\t\t\t[x4,y4],\n")
file_write.write("\t\t\t[x5,height]],\n")
file_write.write("\t\t\twidth=[thk,width2,width3,width4,thk],resolution=6,centered=true);\n")
file_write.write("\t}\n")
file_write.write("\n")
file_write.write("\tmodule lowerPart() {\n")
file_write.write("\t\ttranslate([-10.0,-dia1,0.0])\n")
file_write.write("\t\t%cube([20, thk, 10], center=false);\n")
file_write.write("\t\t%translate([10,-dia1/2.0])\n")
file_write.write("\t\t\tdifference() {\n")
file_write.write("\t\t\t\tcylinder(h=10,d1=dia1,d2=dia1,center=false,$fn=40);\n")
file_write.write("\t\t\t\ttranslate([0,0,-1])\n")
file_write.write("\t\t\t\tcylinder(h=12, d1=dia2, d2=dia2,center=false,$fn=40);\n")
file_write.write("\t\t\t\ttranslate([-2*dia1,-dia1,-1])\n")
file_write.write("\t\t\t\tcube([2*dia1,2*dia1,12], center=false);")
file_write.write("\t\t\t}\n")
file_write.write("\t\ttranslate([-10.0,-thk,0])\n")
file_write.write("\t\tcube([20,thk,10],center=false);\n")
file_write.write("\t}\n")
file_write.write("\n")
file_write.write("\tmodule upperPart() {\n")
file_write.write("\t\ttranslate([-10.0,-dia1,0.0])\n")
file_write.write("\t\tcube([20, thk, 10], center=false);\n")
file_write.write("\t\t%translate([10,-dia1/2,0])\n")
file_write.write("\t\t\tdifference() {\n")
file_write.write("\t\t\t\tcylinder(h=10, d1=dia1, d2=dia1 , center=false, $fn=40);\n")
file_write.write("\t\t\t\ttranslate([0,0,-1])\n")
file_write.write("\t\t\t\tcylinder(h=12, d1=dia2, d2=dia2 , center=false, $fn=40);\n")
file_write.write("\t\t\t\ttranslate([-2*dia1,-dia1,-1])\n")
file_write.write("\t\t\t\tcube([2*dia1,2*dia1,12], center=false);\n")
file_write.write("\t\t\t}\n")
file_write.write("\t\t%translate([-10.0,-thk,0])\n")
file_write.write("\t\tcube([20, thk, 10], center=false);\n")
file_write.write("\t}\n")
file_write.write("}\n")
file_write.close()

print name