import matplotlib.pyplot as plt
dist = []
load = []
lowerbd = []
upperbd = []
#f= open('syntheticGeneratedLoads.txt', 'r')
# f= open('youngsspringneo.txt', 'r')
# f2 = open('youngsspringneoprq7.txt', 'r')
# f3 = open('youngsspringneoprq4.txt', 'r')
# f4 = open('youngsspringneoprq2.5.txt', 'r')
f5 = open('youngsspringneoprq2.txt', 'r')

# for line in f:
# 	step = line.strip("\n").split(",")
# 	dist.append(step[0])
# 	load.append(step[1])
# 	lowerbd.append(step[2])
# 	upperbd.append(step[3])

# dist2 = []
# load2 = []
# for line in f2:
# 	step = line.strip("\n").split(",")
# 	dist2.append(step[0])
# 	load2.append(step[1])

# dist3 = []
# load3 = []
# for line in f3:
# 	step = line.strip("\n").split(",")
# 	dist3.append(step[0])
# 	load3.append(step[1])

# dist4 = []
# load4 = []
# for line in f4:
# 	step = line.strip("\n").split(",")
# 	dist4.append(step[0])
# 	load4.append(step[1])

dist5 = []
load5 = []
for line in f5:
	step = line.strip("\n").split(",")
	dist5.append(step[0])
	load5.append(step[1])

# dist2 = []
# load2 = []
# lowerbd2 = []
# upperbd2 = []
# for line in f2:
# 	step = line.strip("\n").split(",")
# 	dist2.append(step[0])
# 	load2.append(step[1])
# 	lowerbd2.append(step[2])
# 	upperbd2.append(step[3])


# curve = plt.plot(dist2, load2, "ro", label="total")
# curve = plt.plot(dist2, lowerbd2, "--", label="total")
# curve = plt.plot(dist2, upperbd2, "--", label="total")

# plt.plot(dist, load,  "go", label="15k tets")
# plt.plot(dist2, load2, "ro", label ="25k")
# plt.plot(dist3, load3, "bo", label = "35k")
# plt.plot(dist4, load4, "rx", label = "45k")
plt.plot(dist5, load5, "ro", label = "25k")
#curve = plt.plot(dist, lowerbd,  "-", label="total")
#curve = plt.plot(dist, upperbd,  "-", label="total")

plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0.)
plt.show()
