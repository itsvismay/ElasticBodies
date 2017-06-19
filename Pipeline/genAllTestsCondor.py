import sys, os, subprocess

# inputs
#seed = "baseSurrogates/SampleSet_06D_seed01.txt"
seed = "baseSurrogates/SampleSet_02D.txt"
bezierPreped = 'bezierPreped'
bezierMirrorPreped = 'bezierMirrorPreped'
curvePreped = 'curvePreped'
curveMirrorPreped = 'curveMirrorPreped'
ringPreped = 'ringPreped'
origExt = "Orig"
remeshExt = "Remesh"
force = 10000000
experimentDir = "/scratch/cluster/zmisso/ElasticBodies/PipelineTests/"
runF = "/scratch/cluster/zmisso/ElasticBodies/PipelineTests/needToReDo.txt"

data = []
seedFile = open(seed, 'r')
seedLines = seedFile.readlines()
for line in seedLines:
  temp = []
  nums = line.split(',')
  for val in nums:
    temp.append(val.strip())
  data.append(temp)
ind = 0

toRun = []
runFile = open(runF, 'r')
runLines = runFile.readlines()
for line in runLines:
  toRun.append(int(line.strip()))

for row in data:
  print ind, row
  ind = ind + 1

ind = 0

redo = [439, 366, 368, 235, 102, 99, 98, 97, 28, 32, 39, 155, 173, 489, 164, 474, 90, 417, 179, 496, 413, 418, 357, 379, 378, 301, 303, 266, 318, 302, 247, 101, 240, 317, 114, 499, 190, 156, 191, 177, 50, 184, 47, 139, 178, 45, 182, 154, 122, 168, 143, 175, 142, 127, 117, 112, 44, 48, 46, 40, 444, 446, 330, 234]
fails = [238, 406, 270, 347, 308, 219, 189, 83, 138, 85, 11, 5, 460, 407, 384, 344, 276, 238, 225, 188, 165, 110, 36, 30, 494, 382, 355, 320, 283, 222, 185, 172, 118, 66, 107, 451, 419, 373, 331, 270, 212, 134, 193, 176, 129, 68, 440, 408, 388, 342, 252, 218, 153, 115, 119, 67, 51, 416, 363, 358, 304, 246, 210, 195, 94, 197, 120, 17, 395, 405, 354, 322, 209, 147, 148, 162, 132, 78, 49]
#working = [366, 235, 102, 99, 98, 97, 32, 39]

redo2 = [1, 3, 8, 16, 17, 19, 21, 22, 40, 49, 50, 61, 63, 71, 74, 75, 81, 84, 89, 95, 96, 102, 113, 122, 127, 128, 134, 146, 147, 153, 157, 158, 175, 182, 183, 184, 187, 190, 211, 212, 217, 218, 229, 232, 236, 237, 249, 253, 255, 257, 260, 278, 279, 296, 298, 299, 300, 304, 310, 312, 316, 320, 323, 324, 332, 342, 348, 362, 367, 368, 369, 377, 388, 402, 403, 406, 408, 410, 412, 415, 421, 423, 429, 433, 434, 439, 440, 458, 469, 473, 480, 482, 495]
#row = data[563-500]

#for ind in range(1,2):
#for ind in redo:
#for ind in redo2:
for row in data:
#for ind in fails:
#for ind in toRun:
  #row = data[ind]
  print 'Doing for Index', ind
  individualDir = experimentDir + "TwoVarCurveFinalSectionsShort_" + str(ind) + "/"
  
  #subprocess.check_output(['rm', '-rf', individualDir])
  
  mainPipelineCondor = individualDir + 'mainPipelineCondor'
  elasticCondor = individualDir + 'elasticCondor'
  arguements = '/scratch/cluster/zmisso/ElasticBodies/Pipeline/pipeline.py --template /scratch/cluster/zmisso/ElasticBodies/Pipeline/Templates/templateCSpring.py --create ' + individualDir + "optimizeTest.txt --preped " + individualDir + curvePreped + remeshExt + ' --temp ' + individualDir + ' -s'
  
  if not os.path.exists(individualDir):
    os.makedirs(individualDir)
    subprocess.check_output(['python', 'genConfigFile.py', individualDir + "config.txt", individualDir + curvePreped + remeshExt])
    #subprocess.check_output(['python', 'genCurveSpringMirrorData.py', individualDir, str(row[0]), str(row[1]), str(row[2]), str(row[3]), str(row[4]), str(row[5])])
    subprocess.check_output(['python', 'genCurveSpringData.py', individualDir, str(row[0]), "0.6", str(row[1]), "0.6", "0.6"])
    subprocess.check_output(['python', 'generateCondorSubmit.py', '--initialDir', individualDir, '--arguements', arguements, '--file', mainPipelineCondor])
    subprocess.check_output(['python', 'generateCondorSubmit.py', '--initialDir', individualDir, '--arguements', individualDir, '--file', elasticCondor, '--executable', '/scratch/cluster/zmisso/ElasticBodies/Simulation/build/elastic'])
  
    #subprocess.check_output(['python', 'pipeline.py', '--template', 'Templates/templateCSpring.py', '--create', individualDir + 'optimizeTest.txt', '--preped', individualDir + curvePreped + remeshExt, '-s -c'])
  ind = ind + 1

  #print subprocess.check_output(['condor_submit', mainPipelineCondor])
  print str(ind)
  print subprocess.check_output(['condor_submit', elasticCondor])




  #print subprocess.check_output(['python', '/scratch/cluster/zmisso/ElasticBodies/Pipeline/pipeline.py', '--template', '/scratch/cluster/zmisso/ElasticBodies/Pipeline/Templates/templateCSpringMirror.py', '--create', '/scratch/cluster/zmisso/ElasticBodies/Pipeline/optimizeTest.txt', '--sConfig', '/scratch/cluster/zmisso/ElasticBodies/Pipeline/slic3rConfig.ini', '--preped', '/scratch/cluster/zmisso/ElasticBodies/Pipeline/curvePrepedRemesh', '--temp', '/scratch/cluster/zmisso/ElasticBodies/Pipeline', '-s', '-c'])
