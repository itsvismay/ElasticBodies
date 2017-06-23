###############
# LOG METHODS #
###############

import gaobjects

def logFinalResults(population, hallOfFame, experimentDir):
    blah = 5
    # TODO

def logCurrentGeneration(genNumber, experimentDir):
    genFile = open(experimentDir+'generation.txt', 'w')
    genFile.write(str(genNumber))
    genFile.close()

def logid(genNumber, indNumber, indDir, experimentDir):
    individualDir = experimentDir + indDir + "_" + str(genNumber) + "_" + str(indNumber) + "/"
    idFile = open(individualDir+'id.txt', 'w')
    idFile.write(str(indNumber))
    idFile.close()

def logfitness(directory, fit):
    fitFile = open(directory+'fitness.txt', 'w')
    fitFile.write(str(fit))
    fitFile.close()

def printFinalResults(population):
    blah = 5
    # TODO

def printPopulationFitnesses(population, generationNumber):
    for i in range(0, len(population)):
        print population[i].fitness
    print '\n'
