# here goes some stuffs

import sys, getopt, os, subprocess

experimentDir = "/scratch/cluster/zmisso/PipelineTests/"
individualName = "Individual_"
numberOfIndividuals = 10
generationNumber = 0

try:
	opts, args = getopt.getopt(sys.argv[1:], '', ["gen=", "expDir=", "numInd=", "indName="])
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

##################
# NEEDED CLASSES #
##################

class DiscreteVariable:
	value = 0.0
	minValue = 0.0
	maxValue = 0.0
	change = 0.0

	def __init__(self):
		value = 0.0
		# to be implemented

	def __init__(self, val, minVal, maxVal, chg):
		self.value = val
		self.minValue = minVal
		self.maxValue = maxVal
		self.change = chg

class ContinuousVariable:
	value = 0.0
	minValue = 0.0
	maxValue = 0.0

	def __init__(self):
		value = 0.0
		# to be implemented

	def __init__(self, val, minVal, maxVal):
		self.value = val
		self.minValue = minVal
		self.maxValue = maxVal

class Individual:
	continuousVariables = []
	discreteVariables = []

	def __init__(self):
		x = 5
		# to be implemented

	def __init__(self, contVars, discVars):
		self.continuousVariables = contVars
		self.discreteVariables = discVars

	def mutate(self):
		x = 0
		# to be implemented

	def crossover(self, other):
		x = 0
		# to be implemented

#################
# Parse Methods #
#################


######################
# Generation Methods #
######################


########################
# Evolutionary Methods #
########################

def ga(individuals):
	# to be implemented
	start = 5

def mutation(individual):
	number = 5
	# to be implemented

def crossover(individual):
	number = 5
	# to be implemented

def parseIndividuals():
	individuals = []
	# to be implemented
	return individuals

def outputTotalData():
	data = []
	# to be implemented

def outputIterationData():
	data = []
	# to be implemented

###################
# DUMMY FUNCTIONS #
###################

def optimizeForDiscreteValues(individual):
	values = [1, 2, 4, 6, 2]
	# to be implemented

def optimizeForDiscreteFloatValues(individual):
	values = [0.4, 2.8, 0.8, 4.0, 8.0]
	# to be implemented

def optimizeForContinuousValues(individual):
	values = [0.2, 0.7, 0.9, 0.72, 0.3]
	# to be implemented

def optimizeForTwoDiscreteThreeContinuous(individual):
	values = [2, 0.4, 0.39, 0.74, 0.67]
	# to be implemented

##############
# MAIN LOGIC #
##############

toBeImplemented = 0
