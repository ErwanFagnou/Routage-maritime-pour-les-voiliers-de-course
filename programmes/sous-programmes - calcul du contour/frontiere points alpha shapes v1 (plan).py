import random
import matplotlib.pyplot as plt
import numpy as np
from math import *

from scipy.spatial import Delaunay, delaunay_plot_2d

PI = 3.14159265
TWO_PI = 6.28318531

def aff(X, data):
    if data == "POINTS":
        plt.scatter(X[:,0], X[:,1])
    elif data == "LINES":
        plt.plot(X[:,0], X[:,1], 'r')
    plt.show()
    
def dist(A, B):
    return sqrt((A[0]-B[0])**2 + (A[1]-B[1])**2)


def border(X, alpha):
    DT = Delaunay(X[:,:2])

    triangles = DT.simplices
    valid = [False]*len(triangles)
    
    #On enleve les triangles qui ont R>alpha
    for t in range(len(triangles)):
        t1 = triangles[t]
        A, B, C = X[t1]
        a, b, c = dist(B, C), dist(A, C), dist(A, B)
        R = a*b*c / sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c))
        if R <= alpha:
            valid[t] = True

    F = np.empty((0,len(X[0])))
    
    #On regarde quels triangles sont sur un bord, on enregistre les points correspondants
    for t in range(len(triangles)):
        if valid[t]:
            t1 = triangles[t]
            points = X[t1]
            NB = DT.neighbors[t]
            is_border = [(t2==-1 or not valid[t2]) for t2 in NB]
            for i in range(3):
                if is_border[i-1] or is_border[i]:
                    F = np.append(F, [points[i-2]], axis=0)
                
    
    plt.triplot(X[:,0], X[:,1], triangles, color='grey')
    plt.triplot(X[:,0], X[:,1], triangles[valid], color='blue')
    plt.scatter(F[:,0], F[:, 1], color='red')
    plt.show()

    return F


    
n = 2000
alpha = 0.5

X = np.random.normal(0,1,(n,3))
F = border(X, alpha)










