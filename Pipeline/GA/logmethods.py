###############
# LOG METHODS #
###############

import gaobjects

def logFinalResults(population, hallOfFame, experimentDir):
    blah = 5
    # to be implemented

def logPopulationFitnesses(population, generationNumber, experimentDir):
    blah = 5
    # to be implemented

def logPopulationFull(population, generationNumber, experimentDir):
    blah = 5
    # to be implemented

def logIndividualFull(individual, individualNumber, experimentDir):
    blah = 5
    # to be implemented

def logHallOfFame(hallOfFame, experimentDir):
    blah = 3
    # to be implemented

def logCurrentGeneration(genNumber, experimentDir):
    genFile = open(experimentDir+'generation.txt', 'w')
    genFile.write(str(genNumber))
    genFile.close()

def printFinalResults(population):
    blah = 5
    # to be implemented

def printPopulationFitnesses(population, generationNumber):
    for i in range(0, len(population)):
        print population[i].fitness
    print '\n'

def printPopulationFull(population, generationNumber):
    blah = 5
    # to be implemented

def printIndividualFull(individual, individualNumber):
    blah = 5
    # to be implemented
