import matplotlib.pyplot as plt
import sys
t = []
x = []
f = open('redbeamdata.txt')
# filename = "8300000000.000000@0.000001@43257tets@pRa0.6@position.txt"
# f1 = open('../../TestsResults/Damping/'+filename)
filename = sys.argv[1]
f1 = open(filename)
for line in f:
	
	step = line.strip("\n").split(",")
	#t.append(int(step[0]))
	#t.append(((int(step[0]) - 982)/2051.0)-0.08678)#*86/2051)	
	t.append(((int(step[0]) - 982)/2051.0)-0.01)
	x.append(float(step[1]))

t1 = []
x1 = []
for line in f1:
	step = line.strip("\n").split(",")
	t1.append(int(step[0])*1e-5)
	x1.append(float(step[1]))

plt.plot(t[0:], x[0:], "b", label='RealBeamData')
print(t[211])
plt.plot(t1, x1, "g-", label="TestData")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0.)
plt.show()
