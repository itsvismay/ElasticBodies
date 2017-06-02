import sys, os, subprocess

fileName = sys.argv[1] + "optimizeTest.txt"
fileName2 = sys.argv[1] + "points.txt"

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

#minInThk = 1.0
#maxInThk = 5.0
#minInHei = 40.0
#maxInHei = 40.0
#minSections = 4.0
#maxSections = 4.0
#minOutThk = 12.0
#maxOutThk = 12.0
#minWidthIn = 20.0
#maxWidthIn = 20.0

#Data = []
#fr = open(seedData)
#for line in fr.readlines():
#    lineArr = []

#Data = []
#fr = open(seedData)
#for line in fr.readlines():
#    lineArr = []
#    curLine = line.strip().split(',')
#    for i in range(len(curLine)):
#        lineArr.append(curLine[i])
#    Data.append(lineArr)

#Data = [ [ 0.5, 0.5, 0.5, 0.5, 0.5 ] ]
Data = [ float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7]) ]

#print len(Data)

#for i in range(0, len(Data)):
file_write = open(fileName, 'w')
file_write.write(fileName + ".scad "+str(minInThk)+" "+str(maxInThk)+" "+str(minInHei)+" "+str(maxInHei)+" "+str(minSections)+" "+str(maxSections)+" "+str(minOutThk)+" "+str(maxOutThk)+" "+str(minWidthIn)+" "+str(maxWidthIn)+" "+str(Data[0])+" "+str(Data[1])+" "+str(Data[2])+" "+str(Data[3])+" "+str(Data[4])+" "+str(Data[5]))
file_write.close()

file_write = open(fileName2, 'w')
file_write.write(str(Data[0])+" "+str(Data[1])+" "+str(Data[2])+" "+str(Data[3])+" "+str(Data[4])+" "+str(Data[5]))
file_write.close()
#    subprocess.check_output(['python', 'pipeline.py', '--template', 'Templates/templateCSpring.py', '--create', 'optimizeTest.txt', '--sConfig', 'slic3rConfig.ini', '--preped', baseName + str(i), '-s', '-c'])
