import numpy as np
import math
rayleigh1 = 0
rayleigh2 = 0
k = 1
m = 1
c = rayleigh1*m + rayleigh2*k

check = c*c - 4*m*k
root1 = -1*c + math.sqrt(check)/(2*m)
root2 = -1*c - math.sqrt(check)/(2*m)
#x = exp( -c*t/(2m)) *[c1 cos(wt) + c2 sin(wt)], w = sqrt(4mk - c^2)/(2m)

t = 0.00001
for i in range (0, 1000):
    #CASE 1: Overdamped
    if(check > 0):
        c2 = (init_vel + root1*init_pos)/(root1 + root2)
        c1 = init_pos - c2
        x = c1*exp(root1*t*i)+c2*exp(root2*t*i)

    #CASE 2: Critical
    if(check = 0):
        c1 = init_pos
        c2 = init_vel+ (c1*c)/(2*m)
        x = (c1+ c2*t*i)*exp(-1*c*t*i/2/m)

    #Case 3: Underdamped
    if(check < 0):
        w = math.sqrt(4*m*k - c*c)/(2*m)
        c1 = init_pos
        c2 = (init_vel  + c*c1)/w
        x = exp(-1*c*t*i/m/2)*(c1*math.cos(w*t*i) + c2*math.sin(w*t*i))
    print(t*i, x)
