import numpy as np
from scipy.optimize import fmin_cobyla
import sys
import os
import subprocess
import csv
#--------------------------------------
# Input Descriptions
#--------------------------------------
# pathToData               -> the path to the positions files
# pathToAnalyticalConfig   -> the path to the analytical solution
# springType               -> the spring type
#     0 - Curve Spring
#     1 - Mirror Curve Spring
#     2 - Ring Spring
#     3 - Bezier Spring
#     4 - Mirror Bezier Spring
#     5 - Torsion Spring
# simulationConfig         -> the global simulation config file
#--------------------------------------
# Input
#--------------------------------------
pathToConfig = ''
pathToData = ''
pathToAnalyticalConfig = ''
pathToBaseSurrogate = ''
pathToCondorScript = ''
pathToDesignBounds = ''
pathToIterationQueue = ''
baseRandomSampleRate = 0
baseBreadthSampleRate = 0
baseConvergenceRate = 0
baseConvergenceRadius = 0.1
endConvergenceRadius = 0.001
baseBreadthRadius = 0.2
endBreadthRadius = 0.01
springType = -1

# used to tell this script how to execute next command
runType = 0
# 0 -> Run from the beginning and generate the surrogate
# 1 -> Run next iteration of the simulation
# 2 -> Run next iteration of the simulation after x number of results
# 3 -> Run next iteration of the simulation after specific sims complete

try:
	opts, args = getopt.getopt(sys.argv[1:], '', ["simulationConfig=", "runType="])
except: getopt.GetoptError:
	print 'MAJOR ERROR'
	sys.exit(-2)
for opt, arg in opts:
	if opt == 'simulationConfig':
		pathToConfig = arg
	if opt == 'runType':
		runType = int(arg)

simconfig_read = open(pathToConfig, 'r').readlines()
pathToData = simconfig_read[0]
pathToAnalyticalConfig = simconfig_read[1]
pathToBaseSurrogate = simconfig_read[2]
pathToCondorScript = simconfig_read[3]
pathToDesignBounds = simconfig_read[4]
pathToIterationQueue = simconfig_read[5]
baseRandomSampleRate = int(simconfig_read[6])
baseBreadthSampleRate = int(simconfig_read[7])
baseConvergenceRate = int(simconfig_read[8])
baseConvergenceRadius = float(simconfig_read[9])
endConvergenceRadius = float(simconfig_read[10])
baseBreadthRadius = float(simconfig_read[11])
endBreadthRadius = float(simconfig_read[12])
springType = int(simconfig_read[13])

coef = []
#--------------------------------------
# Demo Optimization Function
#--------------------------------------
def sphere(x):
	return np.sum(x**2)
#--------------------------------------
def fitQuadSurface(X, Y):
	n = X.shape[1]
	print 'N ', n
	minNumSamples = (n + 1) * (n + 2) / 2
	Xnew = np.zeros((m, minNumSamples))
	Xnew[:,0] = np.ones(m)
	print Xnew.shape
	#print Xnew
	# Fill in Xnew matrix:
	for i in range(n):
		Xnew[:, i+1] = X[:, i]
		Xnew[:, i+6] = X[:, i]**2

	#print Xnew

	count = 11
	for i in range(n - 1):
		for k in range(i + 1, n):
			Xnew[:,count] = X[:, i] * X[:, k]
			count += 1

	coef = np.linalg.pinv(np.dot(Xnew.T, Xnew))
	coef = np.dot(coef, Xnew.T)
	coef = np.dot(coef, Y)

	return coef
#--------------------------------------
def predict(x, coef):
	x_new = np.zeros(coef.shape)
	x_new[0] = 1.0
	n = len(x)

	for i in range(n):
		x_new[i+1] = x[i]
		x_new[i+6] = x[i]**2
	count = n * 2 + 1

	for i in range(n - 1):
		for k in range(i + 1, n):
			x_new[count] = x[i] * x[k]

	pred = np.sum(coef * x_new)
	return pred
#--------------------------------------
# Cobyla Functions for Finding Min of Surface
#--------------------------------------
def objective(x):
	if np.amax(x) > 1.0:
		return 10000000000.0
	if np.amin(x) < 0.0:
		return 10000000000.0
	return predict(x, coef)
#--------------------------------------
def g0(x):
	return np.amin(x)
#--------------------------------------
def g1(x):
	if np.amax(x) > 1.0:
		return -np.amax(x)
	return np.amax(x)
#--------------------------------------
# Parsing Experiment Data
#--------------------------------------
def parseExperimentData():
	filesInData = os.listdir(pathToData)
	variables = []
	results = []
	completedNames = []
	queueNames = []
	bounds = []
	# parse bounds from design bounds file
	# parse queue names
	for dataFile in filesInData:
		# run file through damp.py
		# add params for test to list of params
		# add results to list of results
		# add names to list of names
	return results, variables, completedNames, queueNames, bounds
#--------------------------------------
# Test Code
#--------------------------------------

fileID = open('OptimalLHS_05D_seed01.txt', 'rb')
reader = csv.reader(fileID, delimiter=',', quoting=csv.QUOTE_NONE)
data = list(reader)
X = np.array(data).astype('double')
m, n = X.shape
print 'X Shape ', X.shape
print X
Y = map(sphere, X)

coef = fitQuadSurface(X, Y)
x = np.ones(n)
y_hat = predict(x, coef)

print y_hat

constraints = [g0, g1]
h_opt = fmin_cobyla(objective, x, constraints, rhoend = 1e-6, maxfun = 1000)

print '\n', h_opt
#print objective(h_opt)
#print g0(h_opt)
#print g1(h_opt)
#--------------------------------------
# Run Pipeline
#--------------------------------------
def runPipeline(x, bounds, runType):
	runPipelineTest(x, bounds, runType)
	# to be implemented
	# generate pipeline config
	# generate condor script
	# run condor script
	return
#--------------------------------------
# Run Pipeline Test
#--------------------------------------
def runPipelineTest(x, bounds, runType):
	# to be implemented
	# generate pipeline config
	# run pipeline
	return
#--------------------------------------
# Build Surrogate
#--------------------------------------
def buildSurrogate(bounds):
	fileID = open(pathToBaseSurrogate, 'rb')
	reader = csv.reader(fileID, delimiter=',', quoting=csv.QUOTE_NONE)
	data = list(reader)

	X = np.array(data).astype('double')
	m, n = X.shape

	for i in range(0, m):
		runPipeline(X[i, :], bounds, 1)

	return
#--------------------------------------
# Check Number Of Iterations
#--------------------------------------
def checkNumberOfIterations(completedNames, queueNames):
	num = 0
	for name in completedNames:
		if name in queueNames:
			num = num + 1
	return num
#--------------------------------------
# Check Last Iterations Done
#--------------------------------------
def checkLastIterationsDone(completedNames, queueNames):
	num = checkNumberOfIterations(completedNames, queueNames)
	if num == len(queueNames):
		return YES
	return NO
#--------------------------------------
# Generate Random Sample
#--------------------------------------
def generateRandomSample(results, variables, bounds):
	size = len(bounds) / 2
	x = np.random.random_sample((size,))
	runPipeline(x, bounds, 3)
	return
#--------------------------------------
# Generate Breadth Sample
#--------------------------------------
def generateBreadthSample(results, variables, bounds):
	size = len(bounds) / 2
	x = np.random.random_sample((size,))
	runPipeline(x, bounds, 4)
	return
#--------------------------------------
# Find Min Of Surface
#--------------------------------------
def findMinOfSurface(results, variables, bounds):
	# TODO :: Replace the below code with last guess
	size = len(bounds) / 2
	initialGuess = np.random.random_sample(size)

	X = np.array(variables).astype('double')
	Y = np.results(variables).astype('double')

	coef = fitQuadSurface(X, Y)
	
	constraints = [g0, g1]
	yHat = fmin_cobyla(objective, initialGuess, constraints, rhoend = 1e-6, maxfun = 1000)

	return yHat
#--------------------------------------
# Calculate Convergence Radius
#--------------------------------------
def calculateConvergenceRadius(startConvRadius, endConvRadius, x, variables, results):
	# to be implemented
	# ask Cem how to do this
	return startConvRadius
#--------------------------------------
# Run Next Iteration
#--------------------------------------
def runNextIteration(results, variables, bounds):
	size = len(bounds) / 2

	for i in range(0, baseRandomSampleRate):
		generateRandomSample(results, variables, bounds)
	for i in range(0, baseBreadthSampleRate):
		generateBreadthSample(results, variables, bounds)

	yHat = findMinOfSurface(results, variables, bounds)

	convergenceRadius = calculateConvergenceRadius(baseConvergenceRadius, endConvergenceRadius, yHat, variables, results)

	positionsToTest = []

	positionsToTest.append(yHat)

	valuesToIgnorePlus = []
	valuesToIgnoreMinus = []

	for i in range(1, baseConvergenceRate):
		for j in range(0, size):
			if not (j in valuesToIngorePlus):
				yHatTempOne = yhat
				yHatTempOne[j] += convergenceRadius * i
				if yHatTempOne[j] > 1.0:
					yHatTempOne[j] = 1.0
					valuesToIgnorePlus.append[j]
				positionsToTest.append(yHatTempOne)

			if not (j in valuesToIgnoreMinus):
				yHatTempTwo = yhat
				yHatTempTwo[j] -= convergenceRadius * i
				if yHatTempTwo[j] < 0.0:
					yHatTempTwo[j] = 0.0
					valuesToIgnoreMinus.append[j]
				positionsToTest.append(yHatTempTwo)

	# TODO :: Handle iteration queue somehow here

	for row in positionsToTest:
		runPipeline(row, bounds, runType)

	return
#--------------------------------------
# Optimization Pipeline
#--------------------------------------
def optimizationPipeline(runType):
	results, variables, completedNames, queueNames, bounds = parseExperimentData()
	if runType == 0:
		buildSurrogate(bounds)
		sys.exit(1)
	if runType == 1:
		# run based on last iteration
		isDone = checkLastIterationDone(completedNames, queueNames)
		if isDone == YES:
			runNextIteration(results, variables, bounds)
		sys.exit(1)
	if runType == 2:
		# run based on number of iterations
		num = checkNumberOfIterations(completedNames, queueNames)
		if num < iterationsToComplete:
			sys.exit(1)
		runNextIteration(results, variables, bounds)
		sys.exit(1)
	if runType == 3:
		# continue random sample
		generateRandomSample(results, variables, bounds)
	if runType == 4:
		# continue breadth sample
		generateBreadthSample(results, variables, bounds)

#--------------------------------------
# Actual Pipeline
#--------------------------------------
#optimizationPipeline(runType)
