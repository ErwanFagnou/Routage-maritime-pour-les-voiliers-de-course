import random
import matplotlib.pyplot as plt
import numpy as np
from math import *

from scipy.spatial import Delaunay
from mpl_toolkits.mplot3d import Axes3D

PI = 3.14159265
TWO_PI = 6.28318531
Rt = 3443.92


def dist(A, B):
    lat1, long1 = A[:2]
    lat2, long2 = B[:2]
    long12 = long2-long1
    p = cos(lat1+lat2)
    q = cos(lat1-lat2)
    s = ((q-p)+(q+p)*cos(long12))/2
    if s>1:                 #erreurs de calculs (ex: s=1.0000000000000002)
        s=1
    elif s<-1:
        s=-1
    return Rt*acos(s)


def to_cartesian(lat, long):
    return np.cos(lat)*np.cos(long), np.sin(lat), np.cos(lat)*np.sin(long)

def plot_sphere(ax):
    ax.set_aspect('equal')

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    r = 0.99
    x = r * np.outer(np.cos(u), np.sin(v))
    y = r * np.outer(np.sin(u), np.sin(v))
    z = r * np.outer(np.ones(np.size(u)), np.cos(v))

    elev = 0.0
    rot = 90 / 180 * PI
    ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='w', linewidth=0, alpha=0.5)

    a = np.array([-np.sin(elev / 180 * np.pi), 0, np.cos(elev / 180 * np.pi)])
    b = np.array([0, 1, 0])
    b = b * np.cos(rot) + np.cross(a, b) * np.sin(rot) + a * np.dot(a, b) * (1 - np.cos(rot))
    ax.plot(np.sin(u),np.cos(u),0,color='k', linestyle = 'dashed')
    horiz_front = np.linspace(0, np.pi, 100)
    ax.plot(np.sin(horiz_front),np.cos(horiz_front),0,color='k')
    vert_front = np.linspace(np.pi / 2, 3 * np.pi / 2, 100)
    ax.plot(a[0] * np.sin(u) + b[0] * np.cos(u), b[1] * np.cos(u), a[2] * np.sin(u) + b[2] * np.cos(u),color='k', linestyle = 'dashed')
    ax.plot(a[0] * np.sin(vert_front) + b[0] * np.cos(vert_front), b[1] * np.cos(vert_front), a[2] * np.sin(vert_front) + b[2] * np.cos(vert_front),color='k')

def aff(X, data):
    if data == "POINTS":
        plt.scatter(X[:,0], X[:,1])
    elif data == "LINES":
        plt.plot(X[:,0], X[:,1], 'r')
    plt.show()
    
def stereo_proj(X):
    proj = np.empty((len(X), 2))
    radius = np.tan(PI/4 + X[:, 0]/2)
    proj[:, 0] = radius * np.cos(X[:, 1])
    proj[:, 1] = radius * np.sin(X[:, 1])
    return proj

def rev_stereo_proj(proj):
    X = np.empty((len(proj), 2))
    radius = np.sqrt(proj[:, 0]**2 + proj[:, 1]**2)
    X[:, 0] = (2*np.arctan(radius)-PI/2)
    X[:, 1] = np.arctan2(proj[:, 1], proj[:, 0])
    return X

def border(X, alpha, pt_out):
    proj = stereo_proj(X)
    
    DT = Delaunay(proj)

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
    bords = np.empty((0,2,len(X[0])))
    #On regarde quels triangles sont sur un bord, on enregistre les points correspondants
    for t in range(len(triangles)):
        if valid[t]:
            t1 = triangles[t]
            points = X[t1]
            NB = DT.neighbors[t]
            is_border = [(t2==-1 or not valid[t2]) for t2 in NB]
            for i in range(3):
                if is_border[i]:
                    bords = np.append(bords, [[points[i-2], points[i-1]]], axis=0)
                if is_border[i-2] or is_border[i-1]:
                    F = np.append(F, [points[i]], axis=0)


    if True:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        plot_sphere(ax)
        
        x, y, z = to_cartesian(X[:,0], X[:,1])
        ax.scatter3D(x, y, z)

        x2, y2, z2 = to_cartesian(F[:,0], F[:,1])
        ax.scatter3D(x2, y2, z2, c='r')
        for b in bords:
            a, b, c = to_cartesian(b[:, 0], b[:, 1])
            plt.plot(a, b, c, color='red')
        plt.show()
        
        plt.triplot(X[:,0], X[:,1], triangles, color='grey')
        plt.triplot(X[:,0], X[:,1], triangles[valid], color='blue')
        plt.scatter(F[:,0], F[:, 1], color='red')
        for b in bords:
            plt.plot(b[:, 0], b[:, 1], color='red')
        plt.show()
    return F

    
n = 2000
alpha = 300

pt_out = [0.5,3]


X = np.random.normal(0,0.3,(n,3))


plt.scatter(X[:,0], X[:,1])
plt.show()


F = border(X, alpha, pt_out)










