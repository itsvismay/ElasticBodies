import matplotlib.pyplot as plt
t = []
x = []
f = open('redbeamdata.txt')

for line in f:
	step = line.strip("\n").split(",")
	t.append(int(step[0])*(86.0/2051))
	x.append(float(step[1]))

plt.plot(t, x, "b-", label='RealBeamData')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0.)
plt.show()
