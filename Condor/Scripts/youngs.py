import matplotlib.pyplot as plt
dist = []
load = []
lowerbd = []
upperbd = []
#f= open('syntheticGeneratedLoads.txt', 'r')
#f = open('youngsnew.txt', 'r')
f = open('youngsspringsvkhalfpoints.txt', 'r')
f2 = open('youngsspringsvk.txt', 'r')

for line in f:
	step = line.strip("\n").split(",")
	dist.append(step[0])
	load.append(step[1])
	lowerbd.append(step[2])
	upperbd.append(step[3])

dist2 = []
load2 = []
lowerbd2 = []
upperbd2 = []
for line in f2:
	step = line.strip("\n").split(",")
	dist2.append(step[0])
	load2.append(step[1])
	lowerbd2.append(step[2])
	upperbd2.append(step[3])


curve = plt.plot(dist2, load2, "ro", label="total")
curve = plt.plot(dist2, lowerbd2, "--", label="total")
curve = plt.plot(dist2, upperbd2, "--", label="total")

curve = plt.plot(dist, load,  "go", label="total")
curve = plt.plot(dist, lowerbd,  "-", label="total")
curve = plt.plot(dist, upperbd,  "-", label="total")

plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0.)
plt.show()
