import matplotlib.pyplot as plt
t=[]
x = []
s = []
k = []
f= open('energy.txt', 'r')
strain= open('senergy.txt', 'r')
kinetic= open('kenergy.txt', 'r')
i=0
maxInt = 70000
for line in f:
	step = line.strip("\n").split(",")
	t.append(step[0])
	x.append(step[1])
	i+=1
	if(i>maxInt):
		break
i=0
for line in strain:
	step = line.strip("\n").split(",")
	s.append(step[1])
	i+=1
	if(i>maxInt):
		break
i=0
for line in kinetic:
	step = line.strip("\n").split(",")
	k.append(step[1])
	i+=1
	if(i>maxInt):
		break

f.close()
total, = plt.plot(t, x, "-", label="total")
strain, = plt.plot(t, s, "ro", label="strain")
kinetic, = plt.plot(t, k, "bo", label = "kinetic")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.)
plt.show()
