import matplotlib.pyplot as plt
t=[]
x = []
f= open('energy.txt', 'r')
i=0
for line in f:
	step = line.strip("\n").split(",")
	t.append(step[0])
	x.append(step[1])
	i+=1
	if(i>10000):
		break
f.close()
plt.plot(t,x, "ro")
plt.show()
