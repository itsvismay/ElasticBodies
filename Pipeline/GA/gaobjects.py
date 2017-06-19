import dummyfitness
import random
import sys

###########################
# DUMMY FITNESS FUNCTIONS #
###########################

import gaobjects

def optimizeForDiscreteValues(individual):
    values = [1, 9, 4, 6, 2]
    discVars = individual.discreteVariables;
    sqerr = 0.0
    sqerr = sqerr + (discVars[0].value - values[0]) * (discVars[0].value - values[0])
    sqerr = sqerr + (discVars[1].value - values[1]) * (discVars[1].value - values[1])
    sqerr = sqerr + (discVars[2].value - values[2]) * (discVars[2].value - values[2])
    sqerr = sqerr + (discVars[3].value - values[3]) * (discVars[3].value - values[3])
    sqerr = sqerr + (discVars[4].value - values[4]) * (discVars[4].value - values[4])
    return sqerr

def optimizeForDiscreteFloatValues(individual):
    values = [0.4, 2.8, 0.8, 4.0, 8.0]
    discVars = individual.discreteVariables;
    sqerr = 0.0
    sqerr = sqerr + (discVars[0].value - values[0]) * (discVars[0].value - values[0])
    sqerr = sqerr + (discVars[1].value - values[1]) * (discVars[1].value - values[1])
    sqerr = sqerr + (discVars[2].value - values[2]) * (discVars[2].value - values[2])
    sqerr = sqerr + (discVars[3].value - values[3]) * (discVars[3].value - values[3])
    sqerr = sqerr + (discVars[4].value - values[4]) * (discVars[4].value - values[4])
    return sqerr

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

def optimizeForTwoDiscreteThreeContinuous(individual):
    values = [2, 0.4, 0.39, 0.74, 0.67]
    contVars = individual.continuousVariables
    discVars = individual.discreteVariables
    sqerr = 0.0
    sqerr = sqerr + (discVars[0].value - values[0]) * (discVars[0].value - values[0])
    sqerr = sqerr + (discVars[1].value - values[1]) * (discVars[1].value - values[1])
    sqerr = sqerr + (contVars[0].value - values[2]) * (contVars[0].value - values[2])
    sqerr = sqerr + (contVars[1].value - values[3]) * (contVars[1].value - values[3])
    sqerr = sqerr + (contVars[2].value - values[4]) * (contVars[2].value - values[4])
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

    def copy(self):
        newVar = DiscreteVariable(self.value, self.minValue, self.maxValue, self.change)
        # newVar.value = self.value
        # newVar.minValue = self.minValue
        # newVar.maxValue = self.maxValue
        # newVar.change = self.change
        return newVar

class ContinuousVariable:
    value = 0.0
    minValue = 0.0
    maxValue = 0.0

    def __init__(self):
        value = 0.0
        minVal = 0.0
        maxVal = 1.0

    def __init__(self, val, minVal, maxVal):
        self.value = val
        self.minValue = minVal
        self.maxValue = maxVal

    def mutate(self):
        if bool(random.getrandbits(1)) == True:
            value = value + random.random() * (maxValue - minValue)
        else:
            value = value - random.random() * (maxValue - minValue)

    def copy(self):
        newVar = ContinuousVariable(self.value, self.minValue, self.maxValue)
        # newVar.value = self.value
        # newVar.minValue = self.minValue
        # newVar.maxValue = self.maxValue
        return newVar

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

    def __repr__(self):
        return "this is an individual"

    def mutate(self, mutationRate):
        x = 0
        newIndividual = self.copy()
        totalVars = len(newIndividual.continuousVariables) + len(newIndividual.discreteVariables)
        for i in range(0, mutationRate):
            varToMutate = random.randint(0, totalVars-1)
            if varToMutate >= len(newIndividual.continuousVariables):
                newIndividual.discreteVariables[varToMutate - len(newIndividual.continuousVariables)].mutate()
            else:
                newIndividual.continuousVariables[varToMutate].mutate()
        return newIndividual

    def crossover(self, other, crossoverRate):
        newIndividual = self.copy()
        totalVars = len(newIndividual.continuousVariables) + len(newIndividual.discreteVariables)
        pointToCrossover = random.randint(0, totalVars-1)
        if pointToCrossover >= len(newIndividual.continuousVariables):
            pos = pointToCrossover - len(newIndividual.continuousVariables)
            newIndividual.discreteVariables[pos] = other.discreteVariables[pos].copy()
        else:
            newIndividual.continuousVariables[pointToCrossover] = other.continuousVariables[pointToCrossover].copy()
        return newIndividual

    def evaluateFitness(self):
        # print 'before fitness eval'
        self.fitness = optimizeForDiscreteValues(self)
        # print 'after fitness eval'
        # fitness = optimizeForDiscreteFloatValues(self)
        # fitness = optimizeForContinuousValues(self)
        # fitness = optimizeForTwoDiscreteThreeContinuous(self)
        if self.fitness == 0.0: # THIS IS TEST CODE REMOVE LATER --- TODO
            print 'FINISHED'
            sys.exit(0)
        return self.fitness

    def copy(self):
        contVars = []
        discVars = []
        for i in range(0, len(self.continuousVariables)):
            contVars.append(self.continuousVariables[i].copy())
        for i in range(0, len(self.discreteVariables)):
            discVars.append(self.discreteVariables[i].copy())
        return Individual(contVars, discVars)

    def getvar(self, index):
        # TODO
        return -1
