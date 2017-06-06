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
