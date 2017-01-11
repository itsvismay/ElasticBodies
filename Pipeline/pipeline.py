# this is the main pipeline script that takes input from the Optimization
# algorithm and transforms the individuals into testable input into the
# simulation.

# input:
# --file + "something.scad"  -> indicates to run with an already generated scad file
# --batch + "filename.txt"   -> indicates to generate scad files for a file containing all of the data for an individual in a single line
# --template + "filename.py" -> indicates which template to use for generating scad files
# --name + "name"            -> indicates specific name for non-batch input file
# --sConfig + "name.ini"     -> indicates the specific slic3r settings to use
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
sConfig = "slic3rConfig.ini"

currentPath = os.path.dirname(os.path.abspath(__file__))

try:
  opts, args = getopt.getopt(sys.argv[1:], 'c', ["file=","batch=","template=","name=","ind=","gen=", "sConfig="])
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
    generation = int(float(arg))
  if opt == "--ind":
    individual = int(float(arg))
  if opt == "-c":
    cleanAll = True
  if opt == "sConfig":
    sConfig = arg

initialScadFiles = []
initialSTLFiles = []
initialGCodeFiles = []
initialLayerSizes = []
layerScadFiles = []
layerSTLFiles = []
layerObjFiles = []
meshedFile = ''
fixedMeshedFile = ''
prepedMesh = ''
forceData = ''

#if isBatch == True:
#  # generate all of the scad files based on the template
#  batchData = open(name)
#  for line in batchData.readlines():
#    curLine = line.strip().split(' ')
#    command = ['python', template, '-a']
#    for obj in curLine:
#      command.append(obj)
#    result = subprocess.check_output(command)
#    initialScadFiles.append(result.strip())
#else:
#  initialScadFiles.append(name)

#for i in range(len(initialScadFiles)):
#  # run openscad -o initialScadFiles[i][:-5]+".stl" initialScadFiles[i]
#
#  print 'openscad -o ' + initialScadFiles[i][:-5] + '.stl ' + initialScadFiles[i]
#  try:
#    result = subprocess.check_output(['openscad', '-o', initialScadFiles[i][:-5] + '.stl', initialScadFiles[i]])
#  except OSError as e:
#    print 'There was a System Error: ', e, '\n'
#  # add generated stl file to list for next step
#  initialSTLFiles.append(initialScadFiles[i][:-5]+".stl")

initialSTLFiles.append("bezspring.stl")

for i in range(len(initialSTLFiles)):
  # run slic3r initialSTLFiles[i] --load sConfig
  print 'slic3r', initialSTLFiles[i], '--load', sConfig
  try:
    result = subprocess.check_output(['slic3r', initialSTLFiles[i], '--load', sConfig])
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
  print 'meshlabserver -i', layerSTLFiles[i], '-o', layerSTLFiles[i][:-4]+'.off'
  try:
    result = subprocess.check_output(['meshlabserver', '-i', layerSTLFiles[i], '-o', layerSTLFiles[i][:-4]+'.off'])
  except OSError as e:
    print 'There was a System Error: ', e, '\n'

  layerObjFiles.append(layerSTLFiles[i][:-4]+'.off')

# combine obj files into one file
for i in range(len(initialGCodeFiles)):
  # run ./3dUnion_bin initialGCodeFiles[i][:-6] initialLayerSizes[i]
  print './3dUnion_bin', initialGCodeFiles[i][:-6], initialLayerSizes[i]
  try:
    temp = 0
    result = subprocess.check_output(['../../libigl/tutorial/build/3dUnion_bin', str(initialGCodeFiles[i][:-6]), str(initialLayerSizes[i])])
  except OSError as e:
    print 'There was a System Error: ', e, '\n'

# fix mesh
meshedFile = "unioned.off"
fixedMeshedFile = "fixedUnion.off"
largestOff = "largest.off"
restOff = "rest.off"
pathToCgal = "../../cgal/Polygon_mesh_processing/examples/Polygon_mesh_processing/"

try:
  print 'split', meshedFile
  result = subprocess.check_output(['./../3dUnion/build/remesh', meshedFile])
except OSError as e:
  print 'There was a System Error', e, '\n'

try:
  print 'remesh', meshedFile
  result = subprocess.check_output(['./'+pathToCgal+'main', largestOff, restOff])
except OSError as e:
  print 'There was a System Error', e, '\n'

# pointless but used for legacy
try:
  print 'cp', 'out.off', fixedMeshedFile
  result = subprocess.check_output(['cp', 'out.off', fixedMeshedFile])
except OSError as e:
  print 'There was a System Error', e, '\n'

# run sim prep to set up the mesh for simulation
force = 10000000;
#prepedMesh = 'prepedMesh.off'
#forceData = 'forcedata.txt'

# temp for beam
prepedMesh = '../shared/beam.off'
forceData = '../shared/beam.txt'
resultMeshes = '../PipelineTests/'

try:
  print './simprep --inputMeshOff', fixedMeshedFile, '--outputM', prepedMesh, 'outputF', forceData, '--forceYAxis', '-1', '--maxForce', force, ''
  result = subprocess.check_output(['./SimPrep/simprep', '--meshOff', fixedMeshedFile, '--outputM', prepedMesh, '--outputF', forceData, '--forceZAxis', '-1', '--maxForce', str(force)])
except OSError as e:
  print 'There was a System Error ', e, '\n'

try:
  print 'cp', fixedMeshedFile, resultMeshes+'IndFixed_'+str(individual)+'.off'
  result = subprocess.check_output(['cp', fixedMeshedFile, resultMeshes+'IndFixed_'+str(individual)+'.off'])
except OSError as e:
  print 'There was a System Error ', e, '\n'

try:
  print 'cp', meshedFile, resultMeshes+'IndUnion_'+str(individual)+'.off'
  result = subprocess.check_output(['cp', meshedFile, resultMeshes+'IndUnion_'+str(individual)+'.off'])
except OSError as e:
  print 'There was a System Error ', e, '\n'

try:
  print 'cp', forceData, resultMeshes+'IndForce_'+str(individual)+'.txt'
  result = subprocess.check_output(['cp', forceData, resultMeshes+'IndForce_'+str(individual)+'.txt'])
except OSError as e:
  print 'There was a System Error ', e, '\n'

# call simulation
#try:
#  print './elastic', '\n'
#  result = subprocess.check_output(['./../elastic'])
#except OSError as e:
#  print 'There was a System Error ', e, '\n'

# lists to clean
# -- initialScadFiles
# -- initialSTLFiles
# -- initialGCodeFiles
# -- layerScadFiles
# -- layerSTLFiles
# -- layerObjFiles
# -- meshedFile
# -- fixedMeshedFile
# -- prepedMesh
# -- forceData

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
  print 'rm', meshedFile
  result = subprocess.check_output(['rm', meshedFile])
  print 'rm', fixedMeshedFile
  result = subprocess.check_output(['rm', fixedMeshedFile])
  #print 'rm', prepedMesh
  #result = subprocess.check_output(['rm', prepedMesh])
  #print 'rm', forceData
  #result = subprocess.check_output(['rm', forceData])
