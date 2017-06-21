#################
# Parse Methods #
#################

import os.path
from gaobjects import DiscreteVariable
from gaobjects import ContinuousVariable
from gaobjects import Individual

# THIS STUFF IS FOR CONDOR CONFIG PARSING AND DATA PARSING

def parseConfig(experimentDir, configName):
    # returns # discrete variables, # continuousVariables, settings
    #
    # settings:
    # [0] - minVal for discrete
    # [1] - maxVal for discrete
    # [2] - minVal for continuous
    # [3] - maxVal for continuous
    #
    config = open(experimentDir+configName, 'r')
    configLines = config.readlines()
    numberOfDiscreteVars = int(configLines[0])
    numberOfContinuousVars = int(configLines[1])
    settings = []
    # TODO -- incorporate change value for discrete vars
    settings.append(float(configLines[2])) # min discrete value
    settings.append(float(configLines[3])) # max discrete value
    settings.append(float(configLines[4])) # min continuous value
    settings.append(float(configLines[5])) # max continuous value
    config.close()

    return numberOfDiscreteVars, numberOfContinuousVars, settings

def parseGenerationNumber(experimentDir):
    if os.path.exists(experimentDir+'generation.txt'):
        genFile = open(experimentDir+'generation.txt', 'r')
        genLines = genFile.readlines()
        generationNumber = int(genLines[0])
        genFile.close()
        return generationNumber
    else:
        return -1

def parseFitness(experimentDir, indDir, genNumber, indNumber):
    individualDir = experimentDir + indDir + "_" + str(genNumber) + "_" + str(indNumber) + "/"
    if os.path.exists(individualDir+'fitness.txt'):
        fitFile = open(individualDir+'fitness.txt', 'r')
        fitLines = fitFile.readlines()
        fitness = float(fitLines[0])
        fitFile.close()
        return fitness
    else:
        return 1000000.0

def parseIndividual(experimentDir, indDir, genNumber, indNumber, settings, numDisc, numCont):
    individualDir = experimentDir + indDir + "_" + str(genNumber) + "_" + str(indNumber) + "/"
    points = []
    discVars = []
    contVars = []
    print 'IN PARSE METHOD'
    popid = -1
    if os.path.exists(individualDir+'points.txt'):
        pointsFile = open(individualDir+'points.txt', 'r')
        for line in pointsFile.readlines():
            points.append(float(line))
        pointsFile.close()
    else:
        print "MAJOR ERROR :: PARSE POPULATION A"
        return
    if os.path.exists(individualDir+"id.txt"):
        idFile = open(individualDir+'id.txt', 'r')
        popid = int(idFile.readlines()[0])
        idFile.close()
    else:
        print "MAJOR ERROR :: PARSE POPULATION B"
        return
    index = 0
    for i in range(0, numDisc): # TODO -- fix change
        discVars.append(DiscreteVariable(points[index], settings[0], settings[1], 1))
        index = index + 1
    for i in range(0, numCont):
        contVars.append(ContinuousVariable(points[index], settings[2], settings[3]))

    if len(points) != len(contVars) + len(discVars):
        print "MAJOR ERROR :: PARSE POPULATION C"
    print len(contVars), 'Number Of Continuous Variables'
    print len(discVars), 'Number Of Discrete Variables'

    return Individual(contVars, discVars, popid);

def parseHallOfFame(experimentDir, hallDir, numHall, settings, numDisc, numCont): # TODO -- clean up
    hallOfFame = []
    hallOfFameDir = experimentDir + hallDir
    if os.path.exists(hallOfFameDir):
        for i in range(0, numHall):
            individualDir = hallOfFameDir + '_' + str(i) + '/'
            points = []
            discVars = []
            contVars = []
            print 'IN PARSE METHOD'
            popid = -1
            fit = 0.0
            if os.path.exists(individualDir+'points.txt'):
                pointsFile = open(individualDir+'points.txt', 'r')
                for line in pointsFile.readlines():
                    points.append(float(line))
                pointsFile.close()
            else:
                # print "MAJOR ERROR :: PARSE POPULATION A"
                return []
            if os.path.exists(individualDir+"fitness.txt"):
                idFile = open(individualDir+'fitness.txt', 'r')
                fit = float(idFile.readlines()[0])
                idFile.close()
            else:
                # print "MAJOR ERROR :: PARSE POPULATION B"
                return []
            index = 0
            for i in range(0, numDisc): # TODO -- fix change
                discVars.append(DiscreteVariable(points[index], settings[0], settings[1], 1))
                index = index + 1
            for i in range(0, numCont):
                contVars.append(ContinuousVariable(points[index], settings[2], settings[3]))

            if len(points) != len(contVars) + len(discVars):
                print "MAJOR ERROR :: PARSE POPULATION C"
            print len(contVars), 'Number Of Continuous Variables'
            print len(discVars), 'Number Of Discrete Variables'
            hallOfFame.append(Individual(contVars, discVars, i))
            hallOfFame[i].fitness = fit
    else:
        os.makedirs(hallOfFameDir)
        for i in range(0, numHall):
            individualDir = hallOfFameDir + '_' + str(i) + '/'
            os.makedirs(individualDir)
    return hallOfFame

def parseSeed(seedPath):
    value = 4
    # TODO
