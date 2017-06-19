#################
# Parse Methods #
#################

import os.path

# THIS STUFF IS FOR CONDOR CONFIG PARSING AND DATA PARSING

def dummyParseConfigOne():
    numberOfDiscreteVars = 5
    numberOfContinuousVars = 0
    settings = [0, 10, 1, 0.0, 1.0]
    return numberOfDiscreteVars, numberOfContinuousVars, settings

def dummyParseConfigTwo():
    numberOfDiscreteVars = 5
    numberOfContinuousVars = 0
    settings = []
    # to be implemented
    return numberOfDiscreteVars, numberOfContinuousVars, settings

def dummyParseConfigThree():
    numberOfDiscreteVars = 5
    numberOfContinuousVars = 0
    settings = []
    # to be implemented
    return numberOfDiscreteVars, numberOfContinuousVars, settings

def dummyParseConfigFour():
    numberOfDiscreteVars = 5
    numberOfContinuousVars = 0
    settings = []
    # to be implemented
    return numberOfDiscreteVars, numberOfContinuousVars, settings

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

def parsePopulation(experimentDir, generationNumber):
    value = generationNumber
    # TODO

def parseSeed(seedPath):
    value = 4
    # TODO

def parseHallOfFame(experimentDir):
    value = 4
    # TODO
