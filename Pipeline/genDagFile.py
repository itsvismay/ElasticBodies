import sys, getopt, os

experimentDir = '/scratch/cluster/zmisso/ElasticBodies/PipelineTests'
individualName = 'Individual_'
gaScript = '/scratch/cluster/zmisso/ElasticBodies/Pipeline/gaScript'
genNumber = -1
numberOfIndividuals = -1

try:
    opts, args = getopt.getopt(sys.argv[1:], '', ["experimentDir=", "individualName=", "gaScript=", "genNumber=", "numIndividuals="])
except getopt.GetoptError:
    print 'INPUT ERROR :: GEN DAG FILE'
for opt, arg in opts:
    if opt == '--experimentDir':
        experimentDir = arg.strip(' \t\n\r')
    elif opt == '--individualName':
        individualName = arg.strip(' \t\n\r')
    elif opt == '--gaScript':
        gaScript = arg.strip(' \t\n\r')
    elif opt == '--genNumber':
        genNumber = int(arg)
    elif opt == '--numIndividuals':
        numberOfIndividuals = int(arg)

filePath = experimentDir + 'dagscript.dag'

dagScript = open(filePath, 'w')
for i in range(0, numberOfIndividuals):
    name = str(genNumber) + '_' + str(i)
    individualDir = experimentDir + individualName + "_" + str(genNumber) + "_" + str(i) + "/"
    # TODO -- make mainPipelineCondor a script param
    dagScript.write("JOB " + name + " " + individualDir + 'mainPipelineCondor' + '\n')

gaName = "GA"
dagScript.write("JOB " + gaName + " " + experimentDir + gaScript + '\n')
for i in range(0, numberOfIndividuals):
    name = str(genNumber) + '_' + str(i)
    dagScript.write("PARENT " + name + " CHILD " + gaName + '\n')
