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
# -c                         -> indicates that this script should delete intermediate files after it runs

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
  batchData = open(name)
  for line in batchData.readlines():
    curLine = line.strip().split(' ')
    command = ['python', template, '-a']
    for obj in curLine:
      command.append(obj)
    result = subprocess.check_output(command)
    initialScadFiles.append(result.strip())
else:
  initialScadFiles.append(name)

for i in range(len(initialScadFiles)):
  # run openscad -o initialScadFiles[i][:-5]+".stl" initialScadFiles[i]
  
  print 'openscad -o ' + initialScadFiles[i][:-5] + '.stl ' + initialScadFiles[i]
  try:
    result = subprocess.check_output(['openscad', '-o', initialScadFiles[i][:-5] + '.stl', initialScadFiles[i]])
  except OSError as e:
    print 'There was a System Error: ', e, '\n'
  # add generated stl file to list for next step
  initialSTLFiles.append(initialScadFiles[i][:-5]+".stl")

for i in range(len(initialSTLFiles)):
  # run slic3r initialSTLFiles[i]
  print 'slic3r', initialSTLFiles[i]
  try:
    result = subprocess.check_output(['slic3r', initialSTLFiles[i]])
    print 'Finished Slic3r'
  except OSError as e:
    print 'There was a System Error: ', e, '\n'
  # add generated gcode file to list for next step
  initialGCodeFiles.append(initialSTLFiles[i][:-4]+".gcode")

for i in range(len(initialGCodeFiles)):
  # run python gcode2layers.py --name initialGCodeFiles[i]
  print 'python gcode2layers.py --name', initialGCodeFiles[i]
  # get number of layers and append it to layers
  numberOfLayers = 0
  try:
    print 'Running gcode2layers.py'
    numberOfLayers = int(subprocess.check_output(['python', 'gcode2layers.py', '--name', initialGCodeFiles[i]])) + 1
  except OSError as e:
    print 'There was a System Error: ', e, '\n'

  # add numberOfLayers to list
  initialLayerSizes.append(numberOfLayers)

  for j in range(numberOfLayers):
    # add generated layer files to list of layers
    layerScadFiles.append(initialGCodeFiles[i][:-6]+"_layer_"+str(j)+".scad")

for i in range(len(layerScadFiles)):
  # run openscad -o layerScadFiles[i][:-5]+".stl" layerScadFiles[i]
  print 'openscad -o', layerScadFiles[i][:-5]+'.stl', layerScadFiles[i]
  try:
    result = subprocess.check_output(['openscad', '-o', layerScadFiles[i][:-5]+'.stl', layerScadFiles[i]])
  except OSError as e:
    print 'There was a System Error: ', e, '\n'

  # add generated file to layerSTLFiles
  layerSTLFiles.append(layerScadFiles[i][:-5] + '.stl')

# convert stl layer files to obj files
for i in range(len(layerSTLFiles)):
  #run meshlabserver -i layerSTLFiles[i] -o layerSTLFiles[i][:-4]+'.obj'
  print 'meshlabserver -i', layerSTLFiles[i], '-o', layerSTLFiles[i][:-4]+'.obj'
  try:
    result = subprocess.check_output(['meshlabserver', '-i', layerSTLFiles[i], '-o', layerSTLFiles[i][:-4]+'.obj'])
  except OSError as e:
    print 'There was a System Error: ', e, '\n'

  layerObjFiles.append(layerSTLFiles[i][:-4]+'.obj')

# combine obj files into one file
for i in range(len(initialGCodeFiles)):
  # run ./3dUnion_bin initialGCodeFiles[i][:-6] initialLayerSizes[i]
  print './3dUnion_bin', initialGCodeFiles[i][:-6], initialLayerSizes[i]
  try:
    temp = 0
    #result = subprocess.check_output(['./../../libigl/tutorial/build/3dUnion_bin', str(initialGCodeFiles[i][:-6]), str(initialLayerSizes[i])])
    #result = subprocess.check_output(['./../../libigl/tutorial/build/3dUnion_bin', str(initialGCodeFiles[i][:-6]), str(3)])
  except OSError as e:
    print 'There was a System Error: ', e, '\n'

  meshedFile = "unioned.obj"

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
    if isBatch == True:
      print 'rm', initialScadFiles[i]
      result = subprocess.check_output(['rm', initialScadFiles[i]])
  for i in range(len(initialSTLFiles)):
    print 'rm', initialSTLFiles[i]
    result = subprocess.check_output(['rm', initialSTLFiles[i]])
  for i in range(len(initialGCodeFiles)):
    print 'rm', initialGCodeFiles[i]
    result = subprocess.check_output(['rm', initialGCodeFiles[i]])
  for i in range(len(layerScadFiles)):
    print 'rm', layerScadFiles[i]
    result = subprocess.check_output(['rm', layerScadFiles[i]])
  for i in range(len(layerSTLFiles)):
    print 'rm', layerSTLFiles[i]
    result = subprocess.check_output(['rm', layerSTLFiles[i]])
  for i in range(len(layerObjFiles)):
    print 'rm', layerObjFiles[i]
    result = subprocess.check_output(['rm', layerObjFiles[i]])
