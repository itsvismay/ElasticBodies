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
    generation = int(arg)
  if opt == "--ind":
    individual = int(arg)
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

# fix mesh
meshedFile = "beam.obj"
fixedMeshedFile = "fixedUnion.obj"
doubleMeshedFile = "doubleFixed.obj"

try:
  print 'python fix_mesh.py', meshedFile
  result = subprocess.check_output(['python', '../../PyMesh/scripts/fix_mesh.py', '--detail', 'low', meshedFile, fixedMeshedFile])
except OSError as e:
  print 'There was a System Error ', e, '\n'

#try:
#  print 'python fix_mesh.py', fixedMeshedFile
#  result = subprocess.check_output(['python', '../../PyMesh/scripts/fix_mesh.py', '--detail', 'low', fixedMeshedFile, doubleMeshedFile])
#except OSError as e:
#  print 'There was a System Error ', e, '\n'

# run sim prep to set up the mesh for simulation
force = 200000;
prepedMesh = 'prepedMesh.obj'
forceData = 'forcedata.txt'

# temp for beam
prepedMesh = '../shared/lowDetailBeam.obj'
forceData = '../shared/lowBeamForce.txt'

try:
  print './simprep --in', fixedMeshedFile, '--out', prepedMesh, '--force', forceData, '--maxForce', force
  result = subprocess.check_output(['./SimPrep/simprep', '--in', fixedMeshedFile, '--out', prepedMesh, '--force', forceData, '--maxForce', str(force)])
except OSError as e:
  print 'There was a System Error ', e, '\n'

# call simulation
try:
  print './elastic', '\n'
  result = subprocess.check_output(['./../elastic'])
except OSError as e:
  print 'There was a System Error ', e, '\n'

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