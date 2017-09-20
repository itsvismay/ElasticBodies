# here goes some stuffs

import sys, getopt, os, subprocess
import random
import gaobjects
import parsemethods
import logmethods

###
# DEFAULT VALUES
###

seed = "baseSurrogates/SampleSet_02D.txt"
experimentDir = "/scratch/cluster/zmisso/ElasticBodies/Pipeline/GA/TestDirectory/"
hallOfFameDir = "HallOfFame/"
individualName = "Individual"
configName = "config.txt"
individualName = "Individual_"
# numberOfIndividuals = 500
numberOfIndividuals = 50
generationNumber = 0
maxGenerations = 100
config = "config.txt"
dummyMode = True
varsToMutate = 3
varsToCross = 2
hallOfFameCount = 4
templateType = 2
runFromConfig = False

###
# Template Type Discriptors
###

# 0  -  bezierSpring
# 1  -  bezierMirrorSpring
# 2  -  curveSpring
# 3  -  curveMirrorSpring
# 4  -  ringSpring
# 5  -  torsionSpring

###
# NEEDED PATHS
###

pathToGenConfigFile = '../genConfigFile.py'
pathToGenerateCondorSubmit = '../generateCondorSubmit.py'
pathToGenCurveSpringData = '../genCurveSpringData.py'
pathToGenCurveMirrorSprintData = 'TODO'
pathToGenBezierSpringData = 'TODO'
pathToGenBezierMirrorSpringData = 'TODO'
pathToGenTorsionSpringData = 'TODO'
pathToGenRingSpringData = 'TODO'

# Parameter Boundary Conditions Paths

pathToBoundaryCurveSpringData = 'TODO'
pathToBoundaryCurveSpringMirrorData = 'TODO'
pathToBoundaryBezierSpringData = 'TODO'
pathToBoundaryBezierSpringMirrorData = 'TODO'
pathToBoundaryTorsionSpringData = 'TODO'
pathToBoundaryRingSpringData = 'TODO'

###
# CLEAN THIS UP LATER
###

bezierPreped = 'bezierPreped'
bezierMirrorPreped = 'bezierMirrorPreped'
curvePreped = 'curvePreped'
curveMirrorPreped = 'curveMirrorPreped'
ringPreped = 'ringPreped'
origExt = "Orig"
remeshExt = "Remesh"

###
# INPUT HANDLING
###

try:
	opts, args = getopt.getopt(sys.argv[1:], 'r', ["gen=", "expDir=", "numInd=", "indName=", "config=", "maxGen=", "mutRate=", "crosRate=", "tempType="])
except getopt.GetoptError:
	print 'Error Bad Input'
	sys.exit(-2)
for opt, arg in opts:
	if opt == "--gen":
		generationNumber = int(arg)
	elif opt == "--indName":
		individualName = arg.strip(' \t\n\r')
	elif opt == "--numInd":
		numberOfIndividuals = int(arg)
	elif opt == "--expDir":
		experimentDir = arg.strip(' \t\n\r')
	elif opt == "--config":
		config = arg.strip(' \t\n\r')
		runFromConfig = True
	elif opt == "--maxGen":
		maxGenerations = int(arg)
	elif opt == "--mutRate":
		varsToMutate = float(arg)
	elif opt == "--crosRate":
		varsToCross = float(arg)
	elif opt == "-r":
		runFromConfig = True
	elif opt == "--tempType":
		templateType = int(arg)

###
# GA LOGIC
###

def generateInitialPopulation(numberOfDiscreteVars, numberOfContinuousVars, settings):
	individuals = []
	for i in range(0, numberOfIndividuals):
		contVars = []
		discVars = []
		for j in range(0, numberOfDiscreteVars):
			minVal = settings[0]
			maxVal = settings[1]
			chg = settings[2]
			val = random.randint(0, (maxVal - minVal) / chg) * chg + minVal
			discVars.append(gaobjects.DiscreteVariable(val, minVal, maxVal, chg))
		for j in range(0, numberOfContinuousVars):
			val = random.random()
			minVal = settings[3]
			maxVal = settings[4]
			contVars.append(gaobjects.ContinuousVariable(val, minVal, maxVal))
		individuals.append(gaobjects.Individual(contVars, discVars, i))
	return individuals

def sortByFitness(population):
	return sorted(population, key=lambda x: x.fitness, reverse=False)

def evaluateFitnessesCondor(population, genNumber):
	# TODO -- Move this to log... maybe
	for individual in population:
		individualDir = experimentDir + individualName + "_" + str(genNumber) + "_" + str(individual.popId) + "/"
		mainPipelineCondor = individualDir + 'mainPipelineCondor'
		# elasticCondor = individualDir + 'elasticCondor' # TODO == remove this
		arguements = '/scratch/cluster/zmisso/ElasticBodies/Pipeline/pipeline.py --template /scratch/cluster/zmisso/ElasticBodies/Pipeline/Templates/templateCSpring.py --create ' + individualDir + "optimizeTest.txt --preped " + individualDir + curvePreped + remeshExt + ' --temp ' + individualDir + ' -s'
		os.makedirs(individualDir)
		subprocess.check_output(['python', pathToGenConfigFile, individualDir + "config.txt", individualDir + curvePreped + remeshExt])
		# TODO -- fill in the correct scripts
		if templateType == 0:
			subprocess.check_output(['python', pathToGenBezierSpringData, individualDir, str(individual.getvar(0)), str(individual.getvar(1)), str(individual.getvar(2)), str(individual.getvar(3)), str(individual.getvar(4))])
		elif templateType == 1:
			subprocess.check_output(['python', pathToGenBezierMirrorSpringData, individualDir, str(individual.getvar(0)), str(individual.getvar(1)), str(individual.getvar(2)), str(individual.getvar(3)), str(individual.getvar(4))])
		elif templateType == 2:
			subprocess.check_output(['python', pathToGenCurveSpringData, individualDir, str(individual.getvar(0)), str(individual.getvar(1)), str(individual.getvar(2)), str(individual.getvar(3)), str(individual.getvar(4))])
		elif templateType == 3:
			subprocess.check_output(['python', pathToGenCurveMirrorSpringData, individualDir, str(individual.getvar(0)), str(individual.getvar(1)), str(individual.getvar(2)), str(individual.getvar(3)), str(individual.getvar(4))])
		elif templateType == 4:
			subprocess.check_output(['python', pathToGenRingSpringData, individualDir, str(individual.getvar(0)), str(individual.getvar(1)), str(individual.getvar(2)), str(individual.getvar(3)), str(individual.getvar(4))])
		elif templateType == 5:
			subprocess.check_output(['python', pathToGenTorsionSpringData, individualDir, str(individual.getvar(0)), str(individual.getvar(1)), str(individual.getvar(2)), str(individual.getvar(3)), str(individual.getvar(4))])

		subprocess.check_output(['python', pathToGenCurveSpringData, individualDir, str(individual.getvar(0)), str(individual.getvar(1)), str(individual.getvar(2)), str(individual.getvar(3)), str(individual.getvar(4))])
		subprocess.check_output(['python', pathToGenerateCondorSubmit, '--initialDir', individualDir, '--arguements', arguements, '--file', mainPipelineCondor])
		logmethods.logid(genNumber, individual.popId, individualName, experimentDir)

		###
		# TEST CODE
		###
		# individual.evaluateFitness()
		# logmethods.logfitness(individualDir, individual.fitness)

	# Create the condor file to rerun this method -- TODO -- replace test code
	subprocess.check_output(['python', '/scratch/cluster/zmisso/ElasticBodies/Pipeline/genGAFile.py']) # TODO -- Fix Params

	###
	# Create a Dag Config File for the Condor Batch Job
	###
	print subprocess.check_output(['python', '/scratch/cluster/zmisso/ElasticBodies/Pipeline/genDagFile.py', '--experimentDir', experimentDir, '--individualName', individualName, '--genNumber', str(genNumber), '--numIndividuals', str(numberOfIndividuals)])
	print subprocess.check_output(['condor_submit_dag', experimentDir + 'dagscript.dag'])

	###
	# Run the Dag Script on Condor
	###
	print subprocess.check_output(['condor_submit', experimentDir + 'dagscript.dag.condor.sub'])

def logHallOfFame(hall):
	# TODO -- Move this to log... maybe
	index = 0
	for individual in hall:
		individualDir = experimentDir + hallOfFameDir + '_' + str(individual.popId) + '/'
		mainPipelineCondor = individualDir + 'mainPipelineCondor'
		arguements = '/scratch/cluster/zmisso/ElasticBodies/Pipeline/pipeline.py --template /scratch/cluster/zmisso/ElasticBodies/Pipeline/Templates/templateCSpring.py --create ' + individualDir + "optimizeTest.txt --preped " + individualDir + curvePreped + remeshExt + ' --temp ' + individualDir + ' -s'
		subprocess.check_output(['python', pathToGenConfigFile, individualDir + "config.txt", individualDir + curvePreped + remeshExt])
		# TODO -- find a good way of organizing the different spring types
		subprocess.check_output(['python', pathToGenCurveSpringData, individualDir, str(individual.getvar(0)), str(individual.getvar(1)), str(individual.getvar(2)), str(individual.getvar(3)), str(individual.getvar(4))])
		subprocess.check_output(['python', pathToGenerateCondorSubmit, '--initialDir', individualDir, '--arguements', arguements, '--file', mainPipelineCondor])
		logmethods.logfitness(individualDir, individual.fitness)
		index = index + 1

def assignFitnesses(population, genNumber):
	for individual in population:
		indNum = individual.popId
		individual.fitness = parsemethods.parseFitness(experimentDir, individualName, genNumber, indNum)

def parsePopulation(genNumber, popSize, settings, numDisc, numCont):
	individuals = []
	for i in range(0, popSize):
		individuals.append(parsemethods.parseIndividual(experimentDir, individualName, genNumber, i, settings, numDisc, numCont))
	return individuals

def checkShouldMutate():
	return bool(random.getrandbits(1))

def reproduce(population, mutationRate):
	newPopulation = []
	for i in range(0, len(population)):
		isMutation = checkShouldMutate()
		if isMutation == True:
			indexToMutate = random.randint(0, max((len(population) / 4) - 1, 0))
			newPopulation.append(population[indexToMutate].mutate(mutationRate, i))
		else:
			other = population[random.randint(0, max((len(population) / 4) - 1, 0))]
			newPopulation.append(population[i].crossover(other, varsToCross, i))
	return newPopulation

def updateHallOfFame(hallOfFame, population):
	newHallOfFame = []
	indexOfFame = 0
	indexOfPop = 0
	while indexOfFame < len(hallOfFame) and indexOfPop < len(population) and len(newHallOfFame) < hallOfFameCount:
		if (hallOfFame[indexOfFame].fitness < population[indexOfPop].fitness):
			newHallOfFame.append(hallOfFame[indexOfFame].copy(indexOfFame + indexOfPop))
			indexOfFame = indexOfFame + 1
		else:
			newHallOfFame.append(population[indexOfPop].copy(indexOfFame + indexOfPop))
			indexOfPop = indexOfPop + 1
	while len(newHallOfFame) < hallOfFameCount:
		newHallOfFame.append(population[indexOfPop].copy(indexOfFame + indexOfPop))
		indexOfPop = indexOfPop + 1
	return newHallOfFame

def updateCurrentGeneration(genNumber):
	genNumber = genNumber + 1
	logmethods.logCurrentGeneration(genNumber, experimentDir)

###
# CONDOR CODEFLOW
###

numberOfDiscreteVars = -1
numberOfContinuousVars = -1

population = []
hallOfFame = []

generationNumber = parsemethods.parseGenerationNumber(experimentDir)
numberOfDiscreteVars, numberOfContinuousVars, settings = parsemethods.parseConfig(experimentDir, configName)

print 'STARTING GENERATION:', generationNumber
print ''

###
# Iterative Genetic Algorithm Implementation
###

if generationNumber == -1:
	population = generateInitialPopulation(numberOfDiscreteVars, numberOfContinuousVars, settings);
	evaluateFitnessesCondor(population, generationNumber+1)
	updateCurrentGeneration(generationNumber)

elif generationNumber < maxGenerations:
	population = parsePopulation(generationNumber, numberOfIndividuals, settings, numberOfDiscreteVars, numberOfContinuousVars)
	assignFitnesses(population, generationNumber)
	population = sortByFitness(population)
	hallOfFame = parsemethods.parseHallOfFame(experimentDir, hallOfFameDir, hallOfFameCount, settings, numberOfDiscreteVars, numberOfContinuousVars)
	hallOfFame = updateHallOfFame(hallOfFame, population)
	logHallOfFame(hallOfFame)

	if generationNumber == maxGenerations:
		logmethods.logFinalResults(population, hallOfFame, experimentDir) # TODO
	else:
		newPopulation = reproduce(population, varsToMutate)
		updateCurrentGeneration(generationNumber)
		evaluateFitnessesCondor(newPopulation, generationNumber+1)

else:
	print 'ALL DONE'

print ''
print 'FINISHED EXECUTING'
