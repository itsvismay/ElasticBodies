# this is the main pipeline script that takes input from the Optimization
# algorithm and transforms the individuals into testable input into the
# simulation.

# input:
# --file + "something.scad"  -> indicates to run with an already generated scad file
# --batch + "filename.txt"   -> indicates to generate scad files for a file containing all of the data for an individual in a single line
# --template + "filename.py" -> indicates which template to use for generating scad files
# --name + "name"            -> indicates specific name for non-batch input file
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
  if opt == "--name":
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
  num = 0
else:
  #initialScadFiles.append(name+"_"+str(generation)+"_"+str(individual))
  initialScadFiles.append(name)

#print 'INITIAL SCAD FILES LEN: ', len(initialScadFiles), '\n'

for i in range(len(initialScadFiles)):
  # run openscad -o initialScadFiles[i][:-5]+".stl" initialScadFiles[i]
  
  print 'openscad -o ' + initialScadFiles[i][:-5] + '.stl ' + initialScadFiles[i]
  try:
    #ls_output = subprocess.check_output(['ls'])
    #print ls_output
    #call(['/bin/base', '-i', '-c' 'ls'])
    #print 'called ls\n'
    #retcode = call(['/bin/base', '-i', '-c', 'openscad -o ' + initialScadFiles[i][:-5] + '.stl ' + initialScadFiles[i]])
    result = subprocess.check_output(['openscad', '-o', initialScadFiles[i][:-5] + '.stl', initialScadFiles[i]])
    #if retcode < 0:
    #  print 'OpenScad Error\n'
  except OSError as e:
    print 'There was a System Error: ', e, '\n'
  # add generated stl file to list for next step
  initialSTLFiles.append(initialScadFiles[i][:-5]+".stl")

#sys.exit(1)

for i in range(len(initialSTLFiles)):
  # run slic3r initialSTLFiles[i]
  print 'slic3r', initialSTLFiles[i]
  try:
    result = subprocess.check_output(['slic3r', initialSTLFiles[i]])
    #retcode = call(['/bin/base', '-i', '-c', 'slic3r ' + initialSTLFiles[i]])
    #if retcode < 0:
    #  print 'Slic3r Error\n'
  except OSError as e:
    print 'There was a System Error: ', e, '\n'
  # add generated gcode file to list for next step
  initialGCodeFiles.append(initialSTLFiles[i][:-4]+".gcode")

#sys.exit(1)

for i in range(len(initialGCodeFiles)):
  # run python gcode2layers.py --name initialGCodeFiles[i]
  print 'python gcode2layers.py --name', initialGCodeFiles[i]
  # get number of layers and append it to layers
  numberOfLayers = 0
  try:
    numberOfLayers = int(subprocess.check_output(['python', 'gcode2layers.py', '--name', initialGCodeFiles[i]]))
    #retcode = call(['/bin/base', '-i', '-c', 'python gcode2layers.py --name ' + initialGCodeFiles[i]])
    #if retcode < 0:
    #  print 'Gcode2Layers Error\n'
    #else:
    #  numberOfLayers = retcode # bad design but whatever
  except OSError as e:
    print 'There was a System Error: ', e, '\n'

  # add numberOfLayers to list
  initialLayerSizes.append(numberOfLayers)

  for j in range(numberOfLayers):
    # add generated layer files to list of layers
    layerScadFiles.append(initialGCodeFiles[i][:-6]+"_layer_"+str(j)+".scad")

#sys.exit(1)

for i in range(len(layerScadFiles)):
  # run openscad -o layerScadFiles[i][:-5]+".stl" layerScadFiles[i]
  print 'openscad -o', layerScadFiles[i][:-5]+'.stl', layerScadFiles[i]
  try:
    result = subprocess.check_output(['openscad', '-o', layerScadFiles[i][:-5]+'.stl', layerScadFiles[i]])
    #retcode = call(['/bin/base', '-i', '-c', 'openscad -o ' + layerScadFiles[i][:-5] + '.stl ' + layerScadFiles[i]])
    #if retcode < 0:
    #  print 'OpenScad Error 2\n'
  except OSError as e:
    print 'There was a System Error: ', e, '\n'

  # add generated file to layerSTLFiles
  layerSTLFiles.append(layerScadFiles[i][:-5] + '.stl')

print 'Done Up to Combine Step'

sys.exit(3)

# convert stl layer files to obj files
# combine obj files into one file
# save it and call simulation

# still a work in progress

# lists to clean
# -- initialScadFiles
# -- initialSTLFiles
# -- initialGCodeFiles
# -- layerScadFiles
# -- layerSTLFiles
# -- layerObjFiles

if cleanAll == True:
  for i in range(len(initialScadFiles)):
    call(['/bin/base', '-i', '-c', 'rm ' + initialScadFiles[i]])
  for i in range(len(initialSTLFiles)):
    call(['/bin/base', '-i', '-c', 'rm ' + initialSTLFiles[i]])
  for i in range(len(initialGCodeFiles)):
    call(['/bin/base', '-i', '-c', 'rm ' + initialGCodeFiles[i]])
  for i in range(len(layerScadFiles)):
    call(['/bin/base', '-i', '-c', 'rm ' + layerScadFiles[i]])
  for i in range(len(layerSTLFiles)):
    call(['/bin/base', '-i', '-c', 'rm ' + layerSTLFiles[i]])
  for i in range(len(layerObjFiles)):
    call(['/bin/base', '-i', '-c', 'rm ' + layerObjFiles[i]])
