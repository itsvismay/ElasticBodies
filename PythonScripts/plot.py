import matplotlib.pyplot as plt
t=[]
x = []
s = []
k = []
g = []
f= open('energy.txt', 'r')
strain= open('senergy.txt', 'r')
kinetic= open('kenergy.txt', 'r')
grav = open('genergy.txt', 'r')
i=0
maxInt = 100000
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
i=0
for line in grav:
	step = line.strip("\n").split(",")
	g.append(step[1])
	i+=1
	if(i>maxInt):
		break

f.close()

total, = plt.plot(t, x, "-", label="total")
strain, = plt.plot(t, s, "ro", label="strain")
kinetic, = plt.plot(t, k, "bo", label = "kinetic")
gravity, = plt.plot(t, g, "go", label = "gravity")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0.)
plt.show()

# xMo=[]
# yMo=[]
# zMo=[]
# t=[]
# m = open("momentum.txt", 'r')
# for line in m:
# 	step = line.strip("\n").split(",")
# 	t.append(step[0])
# 	xMo.append(step[1])
# 	yMo.append(step[2])
# 	zMo.append(step[3])
# m.close()
# plt.plot(t, xMo, "bo", label="xMo")
# plt.plot(t, yMo, "ro", label="yMo")
# plt.plot(t, zMo, "go", label="zMo")
# plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0.)
# plt.show()
