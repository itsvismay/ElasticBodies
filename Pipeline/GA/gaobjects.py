import random
import sys
import numpy

###########################
# DUMMY FITNESS FUNCTIONS #
###########################

def optimizeForContinuousValues(individual):
    values = [0.2, 0.7, 0.9, 0.72, 0.3]
    contVars = individual.continuousVariables;
    sqerr = 0.0
    sqerr = sqerr + (contVars[0].value - values[0]) * (contVars[0].value - values[0])
    sqerr = sqerr + (contVars[1].value - values[1]) * (contVars[1].value - values[1])
    sqerr = sqerr + (contVars[2].value - values[2]) * (contVars[2].value - values[2])
    sqerr = sqerr + (contVars[3].value - values[3]) * (contVars[3].value - values[3])
    sqerr = sqerr + (contVars[4].value - values[4]) * (contVars[4].value - values[4])
    return sqerr

##################
# NEEDED CLASSES #
##################

class DiscreteVariable:
    value = 0.0
    minValue = 0.0
    maxValue = 0.0
    change = 0.0

    def __init__(self):
        self.value = 0.0
        self.minValue = 0.0
        self.maxValue = 0.0
        self.change = 0.0

    def __init__(self, val, minVal, maxVal, chg):
        self.value = val
        self.minValue = minVal
        self.maxValue = maxVal
        self.change = chg

    def mutate(self):
        if bool(random.getrandbits(1)) == True:
            self.value = self.value + self.change
        else:
            self.value = self.value - self.change
        if self.value < self.minVal:
            self.value = self.minVal
        if self.value > self.maxVal:
            self.value = self.maxVal

    def copy(self):
        return DiscreteVariable(self.value, self.minValue, self.maxValue, self.change)

class ContinuousVariable:
    value = 0.0
    minValue = 0.0
    maxValue = 0.0

    def __init__(self):
        self.value = 0.0
        self.minVal = 0.0
        self.maxVal = 1.0

    def __init__(self, val, minVal, maxVal):
        self.value = val
        self.minValue = minVal
        self.maxValue = maxVal

    def mutate(self):
        if bool(random.getrandbits(1)) == True:
            self.value = self.value + random.random() * (self.maxValue - self.minValue) * numpy.random.normal(0.0, 0.1, None)
        else:
            self.value = self.value - random.random() * (self.maxValue - self.minValue) * numpy.random.normal(0.0, 0.1, None)
        if self.value < self.minVal:
            self.value = self.minVal
        if self.value > self.maxVal:
            self.value = self.maxVal

    def copy(self):
        return ContinuousVariable(self.value, self.minValue, self.maxValue)

class Individual:
    continuousVariables = []
    discreteVariables = []
    fitness = 0.0
    popId = -1

    def __init__(self):
        self.fitness = -1.0
        self.continuousVariables = []
        self.discreteVariables = []

    def __init__(self, contVars, discVars, popid):
        self.continuousVariables = contVars
        self.discreteVariables = discVars
        self.fitness = -1.0
        self.popId = popid

    # def __repr__(self):
    #     return "this is an individual"

    def mutate(self, mutationRate, popid):
        x = 0
        newIndividual = self.copy(popid)
        totalVars = len(newIndividual.continuousVariables) + len(newIndividual.discreteVariables)
        for i in range(0, mutationRate):
            varToMutate = random.randint(0, totalVars-1)
            if varToMutate >= len(newIndividual.continuousVariables):
                newIndividual.discreteVariables[varToMutate - len(newIndividual.continuousVariables)].mutate()
            else:
                newIndividual.continuousVariables[varToMutate].mutate()
        return newIndividual

    def crossover(self, other, crossoverRate, popid):
        newIndividual = self.copy(popid)
        totalVars = len(newIndividual.continuousVariables) + len(newIndividual.discreteVariables)
        pointToCrossover = random.randint(0, totalVars-1)
        if pointToCrossover >= len(newIndividual.continuousVariables):
            pos = pointToCrossover - len(newIndividual.continuousVariables)
            newIndividual.discreteVariables[pos] = other.discreteVariables[pos].copy()
        else:
            newIndividual.continuousVariables[pointToCrossover] = other.continuousVariables[pointToCrossover].copy()
        return newIndividual

    # NOTE :: this is test code
    def evaluateFitness(self):
        self.fitness = optimizeForContinuousValues(self)
        return self.fitness

    def copy(self, popid):
        contVars = []
        discVars = []
        for i in range(0, len(self.continuousVariables)):
            contVars.append(self.continuousVariables[i].copy())
        for i in range(0, len(self.discreteVariables)):
            discVars.append(self.discreteVariables[i].copy())
        newIndividual = Individual(contVars, discVars, popid)
        newIndividual.fitness = self.fitness
        return newIndividual

    def getvar(self, index):
        if index < len(self.continuousVariables):
            return self.continuousVariables[index].value
        else:
            return self.discreteVariables[index - len(self.continuousVariables)].value
