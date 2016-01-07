import matplotlib.pyplot as plt
t=[]
x = []
s = []
k = []
f= open('energy.txt', 'r')
strain= open('senergy.txt', 'r')
kinetic= open('kenergy.txt', 'r')
i=0
for line in f:
	step = line.strip("\n").split(",")
	t.append(step[0])
	x.append(step[1])
for line in strain:
	step = line.strip("\n").split(",")
	s.append(step[1])
for line in kinetic:
	step = line.strip("\n").split(",")
	k.append(step[1])

f.close()
plt.plot(t, x, "-", t, s, "ro", t, k, "bo")
plt.show()
