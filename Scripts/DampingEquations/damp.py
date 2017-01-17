import numpy as np
import math
import matplotlib.pyplot as plt


rayleigh1 = 0
rayleigh2 = 0
k = 4
m = 1
c = 4 #rayleigh1*m + rayleigh2*k

init_pos = 1
init_vel = 0
check = c*c - 4*m*k
#x = exp( -c*t/(2m)) *[c1 cos(wt) + c2 sin(wt)], w = sqrt(4mk - c^2)/(2m)
t = 1e-4
sim_seconds = 5

time = []
data = []
for i in range (0, int(sim_seconds/t)):
    #CASE 1: Overdamped
    if(check > 0):
    	root1 = (-1*c + math.sqrt(check))/(2*m)
    	root2 = (-1*c - math.sqrt(check))/(2*m)
        c2 = (init_vel - root1*init_pos)/(root2 - root1)
        c1 = init_pos - c2
        # print(c1, c2, root1, root2)
        x = c1*math.exp(root1*t*i)+c2*math.exp(root2*t*i)

    #CASE 2: Critical
    if(check == 0):
        c1 = init_pos
        c2 = init_vel+ (c1*c)/(2*m)
        print(c1, c2)
        break
        x = (c1+ c2*t*i)*math.exp(-1*c*t*i/2/m)

    #Case 3: Underdamped
    if(check < 0):
        w = math.sqrt(4*m*k - c*c)/(2*m)
        c1 = init_pos
        c2 = (init_vel  + c*c1/(2*m))/w
        x = math.exp(-1*c*t*i/m/2)*(c1*math.cos(w*t*i) + c2*math.sin(w*t*i))
    time.append(t*i)
    data.append(x)

plt.plot(time, data, "-", label="curve")
plt.show()