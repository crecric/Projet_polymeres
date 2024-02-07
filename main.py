import numpy as np
import math
import calcul
import visualisation

tf=100   #total time
t0=0
dt=0.1  #time interval

l=30    #length of polymer (or number of elementary units)

s=int((tf-t0)/dt) #number of steps, size of the array

t=np.linspace(t0,tf,s)
t[0]=t0
y=np.zeros(len(t))
for i in range(t0,len(t)):

    y[i]=calcul.f(t[i]) #length^2

visualisation.graph2D(t,y)