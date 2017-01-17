import matplotlib.pyplot as plt
t = []
x = []
f = open('redbeamdata.txt')
f1 = open('../../TestsResults/Boba/tensile_beam/rawTensile@pRa7.5position.txt')
for line in f:
	
	step = line.strip("\n").split(",")
	#t.append(int(step[0]))
	t.append(((int(step[0]) - 982)/2051.0)-0.08678)#*86/2051)	
	x.append(float(step[1]))

t1 = []
x1 = []
for line in f1:
	step = line.strip("\n").split(",")
	t1.append(int(step[0])*1e-5)
	x1.append(float(step[1]))

plt.plot(t[210:], x[210:], "b", label='RealBeamData')
print(t[211])
plt.plot(t1, x1, "g-", label="TestData")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0.)
plt.show()
