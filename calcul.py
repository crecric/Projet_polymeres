import numpy as np
import random

def polymerize(M,i,j,k):
    val = random.randint(1, 6)
    if val==1 and M[i-1,j,k]==0:
        M[i-1,j,k]=1
        return i-1,j,k
    if val==2 and M[i,j+1,k]==0:
        M[i,j+1,k]=1
        return i,j+1,k
    if val==3 and M[i+1,j,k]==0:
        M[i+1,j,k]=1
        return i+1,j,k
    if val==4 and M[i,j-1,k]==0:
        M[i,j-1,k]=1
        return i,j-1,k
    if val==5 and M[i,j,k-1]==0:
        M[i,j,k-1]=1
        return i,j,k-1
    if val==6 and M[i,j,k+1]==0:
        M[i,j,k+1]=1
        return i,j,k+1
    else:
        return polymerize(M,i,j,k)

