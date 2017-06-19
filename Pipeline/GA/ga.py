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
experimentDir = "TestDirectory/"
individualName = "Individual"
configName = "config.txt"
individualName = "Individual_"
numberOfIndividuals = 500
generationNumber = 0
maxGenerations = 100
config = "config.txt"
dummyMode = True
varsToMutate = 3
varsToCross = 2
hallOfFameCount = 10
runFromConfig = False

###
# NEEDED PATHS
###

pathToGenConfigFile = '../genConfigFile.py'
pathToGenCurveSpringData = '../genCurveSpringData.py'
pathToGenerateCondorSubmit = '../generateCondorSubmit.py'

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
	opts, args = getopt.getopt(sys.argv[1:], 'r', ["gen=", "expDir=", "numInd=", "indName=", "config=", "maxGen=", "mutRate=", "crosRate="])
except getopt.GetoptError:
	print 'Error Bad Input'
	sys.exit(-2)
for opt, arg in opts:
	if opt == "--gen":
		generationNumber = int(arg)
	elif opt == "--indName":
		individualName = arg
	elif opt == "--numInd":
		numberOfIndividuals = int(arg)
	elif opt == "--expDir":
		experimentDir = arg
	elif opt == "--config":
		config = arg
		runFromConfig = True
	elif opt == "--maxGen":
		maxGenerations = int(arg)
	elif opt == "--mutRate":
		varsToMutate = float(arg)
	elif opt == "--crosRate":
		varsToCross = float(arg)
	elif opt == "-r":
		runFromConfig = True

###
# GA LOGIC
###

def generateInitialPopulation(numberOfDiscreteVars, numberOfContinuousVars, settings):
	individuals = []
	for i in range(0, numberOfIndividuals):
		contVars = []
		discVars = []
		for j in range(0, numberOfDiscreteVars):
			val = 0.0 # this should be randomly set not set to zero --- TODO
			minVal = settings[0]
			maxVal = settings[1]
			chg = settings[2]
			discVars.append(gaobjects.DiscreteVariable(val, minVal, maxVal, chg))
		for j in range(0, numberOfContinuousVars):
			val = 0.0 # this should be randomly set not set to zero --- TODO
			minVal = settings[3]
			maxVal = settings[4]
			contVars.append(gaobjects.ContinuousVariable(val, minVal, maxVal))
		individuals.append(gaobjects.Individual(contVars, discVars, i))
	return individuals

def sortByFitness(population):
	# print 'length:', len(population)
	return sorted(population, key=lambda x: x.fitness, reverse=False)

def evaluateFitnesses(population):
	if dummyMode == True:
		# print 'TESTLEN:', len(population)
		for i in range(0, len(population)):
			#print population[i]
			population[i].evaluateFitness()
	else:
		for i in range(0, len(population)):
			population[i].evaluateFitness()
		# more stuffs probably

def evaluateFitnessesCondor(population, directory, indDir, genNumber):
	for individual in population:
		individualDir = directory + indDir + "_" + str(genNumber) + "_" + str(individual.popId) + "/"
		mainPipelineCondor = individualDir + 'mainPipelineCondor'
		# elasticCondor = individualDir + 'elasticCondor' # TODO == remove this
		arguements = '/scratch/cluster/zmisso/ElasticBodies/Pipeline/pipeline.py --template /scratch/cluster/zmisso/ElasticBodies/Pipeline/Templates/templateCSpring.py --create ' + individualDir + "optimizeTest.txt --preped " + individualDir + curvePreped + remeshExt + ' --temp ' + individualDir + ' -s'
		os.makedirs(individualDir)
		subprocess.check_output(['python', pathToGenConfigFile, individualDir + "config.txt", individualDir + curvePreped + remeshExt])
		# TODO -- find a good way of organizing the different spring types
		subprocess.check_output(['python', pathToGenCurveSpringData, individualDir, str(individual.getvar(0)), str(individual.getvar(1)), str(individual.getvar(2)), str(individual.getvar(3)), str(individual.getvar(4))])
		subprocess.check_output(['python', pathToGenerateCondorSubmit, '--initialDir', individualDir, '--arguements', arguements, '--file', mainPipelineCondor])
		#subprocess.check_output(['python', 'generateCondorSubmit.py', '--initialDir', individualDir, '--arguements', individualDir, '--file', elasticCondor, '--executable', '/scratch/cluster/zmisso/ElasticBodies/Simulation/build/elastic'])

		### TODO -- UNCOMMENT THIS FOR FINISHED OPT
		# print subprocess.check_output(['condor_submit', mainPipelineCondor])

def assignFitnesses(population, directory):
	# TODO
	blah = 10

def checkShouldMutate():
	return bool(random.getrandbits(1))

def reproduce(population, mutationRate):
	newPopulation = []
	for i in range(0, len(population)):
		isMutation = checkShouldMutate()
		if isMutation == True:
			indexToMutate = random.randint(0, (len(population) / 4) - 1)
			newPopulation.append(population[indexToMutate].mutate(mutationRate))
		else:
			other = population[random.randint(0, (len(population) / 4) - 1)]
			newPopulation.append(population[i].crossover(other, varsToCross))
	return newPopulation

def updateHallOfFame(hallOfFame, population):
	newHallOfFame = []
	indexOfFame = 0
	indexOfPop = 0
	while indexOfFame < len(hallOfFame) and indexOfPop < len(population) and len(newHallOfFame) < hallOfFameCount:
		if (hallOfFame[indexOfFame].fitness < population[indexOfPop].fitness):
			newHallOfFame.append(hallOfFame[indexOfFame].copy())
			indexOfFame = indexOfFame + 1
		else:
			newHallOfFame.append(population[indexOfPop].copy())
			indexOfPop = indexOfPop + 1
	while len(newHallOfFame) < hallOfFameCount:
		newHallOfFame.append(population[indexOfPop].copy())
	for i in range(0, len(newHallOfFame)):
		newHallOfFame[i].evaluateFitness()
	return newHallOfFame

def updateCurrentGeneration(directory, genNumber):
	genNumber = genNumber + 1
	logmethods.logCurrentGeneration(genNumber, directory)

###
# CONDOR CODEFLOW
###

numberOfDiscreteVars = -1
numberOfContinuousVars = -1

population = []
hallOfFame = []

generationNumber = parsemethods.parseGenerationNumber(experimentDir)

print 'STARTING GENERATION:', generationNumber
print ''

if generationNumber == -1:
	numberOfDiscreteVars, numberOfContinuousVars, settings = parsemethods.parseConfig(experimentDir, configName)
	# print numberOfDiscreteVars, ":: NUM DiscreteVars"
	# print numberOfContinuousVars, ":: NUM ContinuousVars"
	# print settings[0], " :: Min Discrete Value"
	# print settings[1], " :: Max Discrete Value"
	# print settings[2], " :: Min Continuous Value"
	# print settings[3], " :: Max Continuous Value"
	population = generateInitialPopulation(numberOfDiscreteVars, numberOfContinuousVars, settings);
	evaluateFitnessesCondor(population, experimentDir, individualName, generationNumber+1)
	updateCurrentGeneration(experimentDir, generationNumber)

elif generationNumber < maxGenerations:
	population = parsePopulation(experimentDir, generationNumber)
	assignFitnesses(population, experimentDir)
	population = sortByFitness(population)
	hallOfFame = parseHallOfFame(experimentDir)
	hallOfFame = updateHallOfFame(hallOfFame, population)

	if generationNumber == maxGenerations:
		logmethods.logFinalResults(population, hallOfFame, experimentDir)
	else:
		newPopulation = reproduce(population, varsToMutate)
		updateCurrentGeneration(experimentDir, generationNumber)
		evaluateFitnessesCondor(newPopulation, experimentDir, individualName, generationNumber+1)
else:
	print 'ALL DONE'
	# TODO
	# clean up
	# TODO delete generation.txt

print ''
print 'FINISHED EXECUTING'
