import sys, getopt, os

#This program creates a scad file given the input parameters
#
# --name + <fileName>        default: temp.scad
# --thkInPlane + #           default: 2.0
# --widthInPlane + #         default: 40
# --thkOutOfPlane + #        default: 10
# --height                   the height of the spring
# --beamWidth                the width of the middle beam
# --scaleX1 + #              default: 2.0
# --scaleY1 + #              default: 1.0
# --scaleX2 + #              default: 1.85
# --scaleY2 + #              default: 1.0
# --overlap + #              default: 0.5
# --sections + #             default: 3
# -a                         read directly as input

### to be used in full optimizations ###
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
########################################

thkInPlane = 0.5;
widthInPlane = 10;
thkOutOfPlane = 10;
beamWidth = 1.5;
height = 30;
sections = 3;
overlap = 0.5;
scaleX1 = 2.0;
scaleX2 = 2.0;
scaleY1 = 1.0;
scaleY2 = 1.0;
name = "ring_spring.scad"

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
	opts, args = getopt.getopt(sys.argv[1:], 'a' ,["name=", "beamWidth=", "thkInPlane=", "widthInPlane=", "thkOutOfPlane=", "height=", "scaleX1=", "scaleX2=", "scaleY1=", "scaleY2=", "overlap=", "sections=", "minThkInPlane=", "maxThkInPlane", "minWidthInPlane=", "maxWidthInPlane=", "minThkOutOfPlane=", "maxThkOutOfPlane=", "minBeamWidth=", "maxBeamWidth=", "minHeight=", "maxHeight=", "minSections=", "maxSections=", "minOverlap=", "maxOverlap=", "minScaleX1=", "maxScaleX1=", "minScaleX2=", "maxScaleX2=", "minScaleY1=", "maxScaleY1=", "minScaleY2=", "maxScaleY2="])
except getopt.GetoptError:
	print "Error Bad Input"
	sys.exit(-2)
for opt, arg in opts:
	if opt == "--name":
		name = arg
	elif opt == "--thkInPlane":
		thkInPlane = float(arg)
	elif opt == "--widthInPlane":
		widthInPlane = float(arg)
	elif opt == "--thkOutOfPlane":
		thkOutOfPlane = float(arg)
	elif opt == "--height":
		height = float(arg)
	elif opt == "--beamWidth":
		beamWidth = float(arg)
	elif opt == "--scaleX1":
		scaleX1 = float(arg)
	elif opt == "--scaleX2":
		scaleX2 = float(arg)
	elif opt == "--scaleY1":
		scaleY1 = float(arg)
	elif opt == "--scaleY2":
		scaleY2 = float(arg)
	elif opt == "--sections":
		sections = int(float(arg))
	elif opt == "--minThkInPlane":
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
		name = sys.argv[2]
		minThkInPlane = float(sys.argv[3])
		maxThkInPlane = float(sys.argv[4])
		minWidthInPlane = float(sys.argv[5])
		maxWidthInPlane = float(sys.argv[6])
		minThkOutOfPlane = float(sys.argv[7])
		maxThkOutOfPlane = float(sys.argv[8])
		minBeamWidth = float(sys.argv[9])
		maxBeamWidth = float(sys.argv[10])
		minHeight = float(sys.argv[11])
		maxHeight = float(sys.argv[12])
		minSections = float(sys.argv[13])
		maxSections = float(sys.argv[14])
		minOverlap = float(sys.argv[15])
		maxOverlap = float(sys.argv[16])
		minScaleX1 = float(sys.argv[17])
		maxScaleX1 = float(sys.argv[18])
		minScaleX2 = float(sys.argv[19])
		maxScaleX2 = float(sys.argv[20])
		minScaleY1 = float(sys.argv[21])
		maxScaleY1 = float(sys.argv[22])
		minScaleY2 = float(sys.argv[23])
		maxScaleY2 = float(sys.argv[24])
		thkInPlane = calculateThkInPlane(float(sys.argv[25]))
		widthInPlane = calculateWidthInPlane(float(sys.argv[26]))
		thkOutOfPlane = calculateThkOutOfPlane(float(sys.argv[27]))
		beamWidth = calculateBeamWidth(float(sys.argv[28]))
		height = calculateHeight(float(sys.argv[29]))
		sections = calculateSections(float(sys.argv[30]))
		overlap = calculateOverlap(float(sys.argv[31]))
		scaleX1 = calculateScaleX1(float(sys.argv[32]))
		scaleX2 = calculateScaleX2(float(sys.argv[33]))
		scaleY1 = calculateScaleY1(float(sys.argv[34]))
		scaleY2 = calculateScaleY2(float(sys.argv[35]))
	else:
		print 'Error UnRecognized Command'
		exit(3)

def calculateThkInPlane(v):
	return (maxThkInPlane - minThkInPlane) * v + minThkInPlane

def calculateWidthInPlane(v):
	return (maxWidthInPlane - minWidthInPlane) * v + minWidthInPlane

def calculateThkOutOfPlane(v):
	return (maxThkOutOfPlane - minThkOutOfPlane) * v + minThkOutOfPlane

def calculateBeamWidth(v):
	return (maxBeamWidth - minBeamWidth) * v + minBeamWidth

def calculateHeight(v):
	return (maxHeight - minHeight) * v + minHeight

def calculateSections(v):
	return math.floor( (maxSections - minSections) * v + minSections )

def calculateOverlap(v):
	return (maxOverlap - minOverlap) * v + minOverlap

def calculateScaleX1(v):
	return (maxScaleX1 - minScaleX1) * v + minScaleX1

def calculateScaleX2(v):
	return (maxScaleX2 - minScaleX2) * v + minScaleX2

def calculateScaleY1(v):
	return (maxScaleY1 - minScaleY1) * v + minScaleY1

def calculateScaleY2(v):
	return (maxScaleY2 - minScaleY2) * v + minScaleY2

def writeSectionData(file_write, sectionNum):
	file_write.write("// "+str(sectionNum+1)+".SECTION\n")
	if sectionNum == 0:
		file_write.write("translate([0.0, -beamWidth + 0.5, 0.0]) cube([widthInPlane,beamWidth,thkOutOfPlane], center=false);")
	if sectionNum != 0:
		file_write.write("translate([0,"+str(sectionNum)+"*(diaExt-overlap),0]){\n")
	file_write.write("cube([lengthOfCube,beamWidth / 2.0,thkOutOfPlane], center=false);\n")
	if sectionNum == int(sections)-1:
		file_write.write("translate([0.0, diaExt - beamWidth / 4.0,0.0]) cube([widthInPlane,beamWidth,thkOutOfPlane], center=false);\n")
	file_write.write("translate([lengthOfCube,0,0]) scaled_half_circle(thkOutOfPlane, diaExt, diaInt);\n")
	file_write.write("translate([0,diaExt- beamWidth / 2.0,0]) cube([lengthOfCube,beamWidth / 2.0,thkOutOfPlane], center=false);\n")
	if sectionNum != 0:
		file_write.write("}\n")
	file_write.write("\n")

#code that makes basic scad file
file_write = open(name, 'w')
file_write.write("\n")
file_write.write("thkInPlane = %3.8f;\n" % thkInPlane)
file_write.write("widthInPlane = %3.8f;\n" % widthInPlane)
file_write.write("thkOutOfPlane = %3.8f;\n" % thkOutOfPlane)
file_write.write("beamWidth = %3.8f;\n" % beamWidth)
file_write.write("\n")
file_write.write("$fn=40;\n")
file_write.write("height = %3.8f;\n" % height)
file_write.write("sections = %d;\n" % sections)
file_write.write("diaExt = height / sections;\n")
file_write.write("scaleX1 = %3.8f;\n" % scaleX1)
file_write.write("scaleY1 = %3.8f;\n" % scaleY1)
file_write.write("scaleX2 = %3.8f;\n" % scaleX2)
file_write.write("scaleY2 = %3.8f;\n" % scaleY2)
file_write.write("overlap = %3.8f;\n" % overlap)
file_write.write("\n")
file_write.write("diaInt = diaExt - 2*thkInPlane;\n")
file_write.write("\n")
file_write.write("lengthOfCube = widthInPlane/2; echo(lengthOfCube);\n")
file_write.write("\n")
file_write.write("module scaled_half_circle(thkOut, dia1, dia2) {\n")
file_write.write("translate([0,dia1/2,0])\n")
file_write.write("intersection(){\n")
file_write.write("difference() {\n")
file_write.write("\tscale([scaleX1, scaleY1, 1]) cylinder(h=thkOut, d1=dia1, d2=dia1, center=false);\n")
file_write.write("\tscale([scaleX2, scaleY2, 1]) translate([0,0,-2.5]) cylinder(h=thkOut+5, d1=dia2, d2=dia2, center=false);\n")
file_write.write("\t}\n")
file_write.write("\ttranslate([dia1/2,0,thkOut/2]) cube([dia1,dia1,thkOut*2], center=true);\n")
file_write.write("\t}\n")
file_write.write("}\n")
file_write.write("\n")
file_write.write("\n")
file_write.write("union() {\n")
for i in range(0,sections):
	writeSectionData(file_write, i)
file_write.write("\n")
file_write.write("mirror([1,0,0]){\n")
for i in range(0,sections):
	writeSectionData(file_write, i)
file_write.write("}\n")
file_write.write("}\n")
file_write.close()

print name
