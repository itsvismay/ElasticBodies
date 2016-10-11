import sys, getopt, os

#This program creates a scad file given the input parameters
#
# -f - full process
# --thkInPlane + #           default: 2.0
# --widthInPlane + #         default: 40
# --thkOutOfPlane + #        default: 10
# --n + #                    default: 40
# --diaExt + #               default: 30
# --scaleX1 + #              default: 2.0
# --scaleY1 + #              default: 1.0
# --scaleX2 + #              default: 1.85
# --scaleY2 + #              default: 1.0
# --overlap + #              default: 0.5
# --sections + #             default: 3
# --name + <fileName>        default: temp.scad

vThkInPlane = 2.0
vWidthInPlane = 40
vThkOutOfPlane = 10
vN = 40
vDiaExt = 30
vScaleX1 = 2.0
vScaleY1 = 1.0
vScaleX2 = 1.85
vScaleY2 = 1.0
vOverlap = 0.5
vSections = 3
vName = 'temp.scad'
vMakeFull = False

#test code
#print 'Number of Arguements:', len(sys.argv), 'arguements.'
#print 'Arguement List:', str(sys.argv)

#parse arguements
try:
	opts, args = getopt.getopt(sys.argv[1:],'f' ,["thkInPlane=","widthInPlane=","thkOutOfPlane=","n=","diaExt=","scaleX1=","scaleY1=","scaleX2=","scaleY2=","overlap=","name=","sections="])
except getopt.GetoptError:
	print "Error Bad Input"
	sys.exit(-2)
for opt, arg in opts:
	if opt == "-f":
		print 'Making Full'
		vMakeFull = True
	elif opt == "--thkInPlane":
		vThkInPlane = arg
	elif opt == "--widthInPlane":
		vWidthInPlane = arg
	elif opt == "--thkOutOfPlane":
		vThkOutOfPlane = arg
	elif opt == "--n":
		vN = arg
	elif opt == "--diaExt":
		vDiaExt = arg
	elif opt == "--scaleX1":
		vScaleX1 = arg
	elif opt == "--scaleY1":
		vScaleY1 = arg
	elif opt == "--scaleX2":
		vScaleX2 = arg
	elif opt == "--scaleY2":
		vScaleY2 = arg
	elif opt == "--overlap":
		vOverlap = arg
	elif opt == "--name":
			vName = arg
	elif opt == "--sections":
			vSections = int(arg)

def writeSectionData(file_write, sectionNum):
	file_write.write("// "+str(sectionNum+1)+".SECTION\n")
	if sectionNum != 0:
		file_write.write("translate([0,"+str(sectionNum)+"*(diaExt-overlap),0]){\n")
	file_write.write("cube([lengthOfCube,2,thkOutOfPlane], center=false);\n")
	file_write.write("translate([lengthOfCube,0,0]) scaled_half_circle(thkOutOfPlane, diaExt, diaInt, scaleX1, scaleY1, scaleX2, scaleY2);\n")
	file_write.write("translate([0,diaExt-2,0]) cube([lengthOfCube,2,thkOutOfPlane], center=false);\n")
	if sectionNum != 0:
		file_write.write("}\n")
	file_write.write("\n")

#code that makes basic scad file
file_write = open(vName, 'w')
file_write.write("///////////////////////////////////////////////////////\n")
file_write.write("// PARAMETERS\n")
file_write.write("thkInPlane = "+str(vThkInPlane)+"; //mm Thickness of the circular section before scaling\n")
file_write.write("widthInPlane = 40; //width of spring in mm\n")
file_write.write("thkOutOfPlane = 10; //mm\n")
file_write.write("\n")
file_write.write("$fn=40;\n")
file_write.write("diaExt = 30;\n")
file_write.write("scaleX1 = 2.0;\n")
file_write.write("scaleY1 = 1.0;\n")
file_write.write("scaleX2 = 1.85;\n")
file_write.write("scaleY2 = 1.0;\n")
file_write.write("overlap = 0.5;\n")
file_write.write("\n")
file_write.write("diaInt = diaExt - 2*thkInPlane;\n")
file_write.write("\n")
file_write.write("//To calculate the length of the bottom flat piece (cube), we need to use \"scaleX\" information:\n")
file_write.write("lengthOfCube = widthInPlane/2 - scaleX1*diaExt/4; echo(lengthOfCube);\n")
file_write.write("\n")
file_write.write("module scaled_half_circle(thkOut, dia1, dia2, scale_x1, scale_y1, scale_x2, scale_y2) {\n")
file_write.write("translate([0,dia1/2,0])\n")
file_write.write("intersection(){\n")
file_write.write("difference() {\n")
file_write.write("\tscale([scale_x1, scale_y1, 1]) cylinder(h=thkOut, d1=dia1, d2=dia1, center=false);\n")
file_write.write("\tscale([scale_x2, scale_y2, 1]) translate([0,0,-2.5]) cylinder(h=thkOut+5, d1=dia2, d2=dia2, center=false);\n")
file_write.write("\t}\n")
file_write.write("\ttranslate([dia1/2,0,thkOut/2]) cube([dia1,dia1,thkOut*2], center=true);\n")
file_write.write("\t}\n")
file_write.write("}\n")
file_write.write("\n")
file_write.write("\n")
for i in range(0,vSections):
	writeSectionData(file_write, i)
file_write.write("// MIRROR\n")
file_write.write("\n")
file_write.write("mirror([1,0,0]){\n")
for i in range(0,vSections):
	writeSectionData(file_write, i)
file_write.write("}\n")
file_write.close()

print "Created Initial Scad File:", str(vName)

import subprocess
from subprocess import call

currentPath = os.path.dirname(os.path.abspath(__file__))

if vMakeFull == True:
	# run openSCAD to create the stl file
	try:
		print 'Calling:', 'openscad -o', vName[:-5]+'.stl', vName, '\n'
		#retcode = call("openscad -o " + currentPath+"/"+vName[:-5]+".stl" + currentPath+"/"+vName, shell=False)
		#retcode = call("openscad", shell=True, executable='/bin/bash')
		retcode = call(['/bin/bash', '-i', '-c', 'openscad -o ' + vName[:-5]+'.stl ' + vName])
		if retcode < 0:
			print 'Script Error\n'
		else:
			print 'Finished Running OpenSCAD\n'
	except OSError as e:
		print 'There Was A System Error:', e, '\n'

	# run slic3r to generate the gcode
	try:
		print 'Calling:', 'ls '+vName[:-5]+'.stl', '\n'
		call(['/bin/bash', '-i', '-c', 'ls'])
		print 'Finished Call\n'
		if retcode < 0:
			print 'Script Error\n'
		else:
			print 'Finished Running Slic3r\n'
	except OSError as e:
		print 'There Was A System Error:', e, '\n'
	finally:
		print 'YEAEHAHSDHF\n'

	# run Cem's script to change gcode into scad layers

print 'ALL DONE'

#temp openscad -o vName[:-5]+".stl" vName
#temp slic3r vName[:-5]+".stl" #config options
#temp python gcode2openscad.py
