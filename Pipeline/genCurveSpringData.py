import sys, os, subprocess

baseName = 'Individual'
seedData = 'OptimalLHS_05D_seed01.txt'
fileName = 'optimizeTest.txt'

minInThk = 0.4
maxInThk = 4.0
minInHei = 12.0
maxInHei = 38.0
minSections = 4.0
maxSections = 10.0
minOutThk = 9.0
maxOutThk = 12.0
minWidthIn = 9.0
maxWidthIn = 12.0

Data = []
fr = open(seedData)
for line in fr.readlines():
    lineArr = []
    curLine = line.strip().split(',')
    for i in range(len(curLine)):
        lineArr.append(curLine[i])
    Data.append(lineArr)

#print len(Data)

#for i in range(0, len(Data)):
for i in range(0, len(Data)):
    file_write = open(fileName, 'w')
    file_write.write(fileName + ".scad "+str(minInThk)+" "+str(maxInThk)+" "+str(minInHei)+" "+str(maxInHei)+" "+str(minSections)+" "+str(maxSections)+" "+str(minOutThk)+" "+str(maxOutThk)+" "+str(minWidthIn)+" "+str(maxWidthIn)+" "+str(Data[i][0])+" "+str(Data[i][1])+" "+str(Data[i][2])+" "+str(Data[i][3])+" "+str(Data[i][4]))
    file_write.close()

    subprocess.check_output(['python', 'pipeline.py', '--template', 'Templates/templateCSpring.py', '--create', 'optimizeTest.txt', '--sConfig', 'slic3rConfig.ini', '--preped', baseName + str(i), '-s', '-c'])
