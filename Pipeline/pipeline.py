# this is the main pipeline script that takes input from the Optimization
# algorithm and transforms the individuals into testable input into the
# simulation.

# input:
# --file + "something.scad"  -> indicates to run with an already generated scad file
# --batch + "filename.txt"   -> indicates to generate scad files for a file containing all of the data for an individual in a single line
# --template + "filename.py" -> indicates which template to use for generating scad files
# --name + "name"            -> indicates specific name for template
# --ind + #                  -> indicates individual number
# --gen + #                  -> indicates generation number
# -c                         -> indicates that this script should delete intermediate files as it runs

# if both file and batch are defined the results are undefined
# name, ind, gen are completely optional
# ind just defines the index of the first individual

import sys, getopt, os, subprocess
from subprocess import call

isBatch = False
isInd = False
cleanAll = False
template = ""
name = ""
generation = 0
individual = 0

currentPath = os.path.dirname(os.path.abspath(__file__))

try:
  opts, args = getopt.getopt(sys.argv[1:], 'c', ["file=","batch=","template=","name=","ind=","gen="])
except getopt.GetoptError:
  print 'Error Bad Input'
  sys.exit(-2)
for opt, arg in opts:
  if opt == "--file":
    isBatch = False
    isFile = True
    name = arg
  if opt == "--batch":
    isBatch = True
    isFile = False
    name = arg
  if opt == "--template":
    template = arg
  if opt == "--name"
    name = arg
  if opt == "--gen":
    generation = int(arg)
  if opt == "--ind":
    individual = int(arg)
  if opt == "-c":
    cleanAll = True

initialScadFiles = []
initialSTLFiles = []
initialGCodeFiles = []
initialLayerSizes = []
layerScadFiles = []
layerSTLFiles = []
layerObjFiles = []
meshedFile = ''

if isBatch == True:
  # generate all of the scad files based on the template
else:
  initialScadFiles.append(name+"_"+str(generation)+"_"+str(individual))

for i in range(len(initialScadFiles)):
  # run openscad -o initialScadFiles[i][:-5]+".stl" initialScadFiles[i]

  try:
    retcode = call(['/bin/base', '-i', '-c', 'openscad -o ' + initialScadFiles[i][:-5] + '.stl ' + initialScadFiles[i]])
    if retcode < 0:
      print 'OpenScad Error\n'
  except OSError as e:
    print 'There was a System Error: ', e, '\n'

  initialSTLFiles.append(initialScadFiles[i][:-5]+".stl")

for i in range(len(initialSTLFiles)):
  # run slic3r initialSTLFiles[i]
  # to be implemented
  initialGCodeFiles.append(initialSTLFiles[i][:-4]+".gcode")

for i in range(len(initialGCodeFiles)):
  # run gcode2layers.py --name initialGCodeFiles[i]
  # get number of layers and append it to layers
  # to be implemented

# change this
for i in range(len(initialGCodeFiles)):
  for j in range(initialLayerSizes[i]):
    # run openscad -o initialGCodeFiles[i][:-5]+"_layer"+str(j)+".stl" initialGCodeFiles[i][:-5]+"_layer"+str(j)+".scad"

# still a work in progress

if cleanAll == True:
  # to be implemented
