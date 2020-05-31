import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from math import *
import time

PI = np.pi
TWO_PI = 2*PI

def dist(A, B):
    lat1 = A[0]
    lat2 = B[0]
    long12 = B[1] - A[1]
    return 60*180/PI*acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(long12))

def angle_rot(A, B):
    lat1 = A[0]
    lat2 = B[0]
    long12 = B[1] - A[1]
    
    return (PI/2-atan2(cos(lat2)*sin(long12),(cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(long12))))%TWO_PI

def is_between_angles(a, a_min, a_max):
    if a_min <= a_max:
        return a_min < a and a < a_max
    else:
        return a_min < a or a < a_max

##def do_segments_intersect(A, B, C, D):
##    a1 = angle_rot(A, B)
##    a2 = angle_rot(C, D)
##    if (a1-a2)%TWO_PI <= PI:
##        return is_between_angles(a1, angle_rot(A, D), angle_rot(A, C)) and is_between_angles(a2, angle_rot(C, A), angle_rot(C, B)) and is_between_angles(angle_rot(B, A) , angle_rot(B, C), angle_rot(B, D)) and is_between_angles(angle_rot(D, C), angle_rot(D, B), angle_rot(D, A))
##    else:
##        return is_between_angles(a1, angle_rot(A, C), angle_rot(A, D)) and is_between_angles(a2, angle_rot(C, B), angle_rot(C, A)) and is_between_angles(angle_rot(B, A) , angle_rot(B, D), angle_rot(B, C)) and is_between_angles(angle_rot(D, C), angle_rot(D, A), angle_rot(D, B))


def do_segments_intersect(A, B, C, D):
    Ac = np.array(to_cartesian(A[0], A[1]))
    Bc = np.array(to_cartesian(B[0], B[1]))
    Cc = np.array(to_cartesian(C[0], C[1]))
    Dc = np.array(to_cartesian(D[0], D[1]))

    V2 = np.cross(Cc, Dc)
    if np.dot(V2, Ac)*np.dot(V2, Bc)<0: #ie A et B de chaque coté du plan (si =0 on ne compte pas l'intersection)
        V1 = np.cross(Ac, Bc)
        if np.dot(V1, Cc)*np.dot(V1, Dc)<0: #idem pour C et D
            V3 = np.cross(Ac, Cc)
            if np.dot(V3, Bc)*np.dot(V3, Dc)>0: #ie si B et D appartiennent au meme hemisphere (A et C dedans car le définissent)
                return True
    return False
    
    

def border(X, pt_out, min_dist, max_dist, show=False):
    best_dist = dist(X[0], pt_out)
    first_pt = X[0]
    for i in range(1, n):
        d = dist(X[i], pt_out)
        if d < best_dist:
            best_dist = d
            first_pt = X[i]

    F = [first_pt]
    prev_angle = angle_rot(first_pt, pt_out) % TWO_PI
    while len(F)==1 or ((F[-1]!=F[0]).any() and len(F)==2) or (((F[-1]!=F[1]).any() or (F[-2]!=F[0]).any()) and len(F)>=3):
        best_angle = TWO_PI + 1
        next_pt = F[-1]
        for p in X:
            if (p != F[-1]).all() and dist(F[-1], p) <= max_dist:
                angle = (angle_rot(F[-1], p) - prev_angle) % TWO_PI
                if angle == 0:
                    angle = TWO_PI
                if angle < best_angle:
                    intersect = do_segments_intersect(F[-1], p, pt_out, first_pt)
                    for i in range(len(F)-1):
                        intersect = intersect or do_segments_intersect(F[-1], p, F[i], F[i+1])
                    if not intersect:
                        best_angle = angle
                        next_pt = p
        prev_angle = angle_rot(next_pt, F[-1])
        F.append(next_pt)
        if len(F)>100:
            print("erreur")
            break
    if show or len(F)>100:
        F2 = np.array(F)
        plt.scatter(X[:,0], X[:,1])
        plt.scatter(F2[:,0], F2[:,1], c='r')
        plt.scatter([pt_out[0]], [pt_out[1]], c='g')
        
        plt.plot([pt_out[0], first_pt[0]], [pt_out[1], first_pt[1]], 'g')
        plt.plot(F2[:,0], F2[:,1], 'r')
        plt.show()

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        plot_sphere(ax)
        
        x, y, z = to_cartesian(X[:,0], X[:,1])
        ax.scatter3D(x, y, z)

        x2, y2, z2 = to_cartesian(F2[:,0], F2[:,1])
        ax.scatter3D(x2, y2, z2, c='r')
        ax.plot3D(x2, y2, z2, 'r')

        ax.scatter3D([cos(pt_out[0])*cos(pt_out[1])], [sin(pt_out[0])], [cos(pt_out[0])*sin(pt_out[1])], c='g')
        ax.plot3D([cos(pt_out[0])*cos(pt_out[1]), cos(first_pt[0])*cos(first_pt[1])], [sin(pt_out[0]), sin(first_pt[0])], [cos(pt_out[0])*sin(pt_out[1]), cos(first_pt[0])*sin(first_pt[1])], c='g')
        plt.show()
        
    return F

def to_cartesian(lat, long):
    return np.cos(lat)*np.cos(long), np.sin(lat), np.cos(lat)*np.sin(long)

def plot_sphere(ax):
    ax.set_aspect('equal')

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    x = 1 * np.outer(np.cos(u), np.sin(v))
    y = 1 * np.outer(np.sin(u), np.sin(v))
    z = 1 * np.outer(np.ones(np.size(u)), np.cos(v))

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


n = 1000
min_dist = 0.0
max_dist = 2000
pt_out = [PI,0]


#X = np.random.random((n,2))
X = np.random.normal(0,1,(n,2))
#X = (2*X-1)/((0.5+np.abs(2*X-1)**3))
X[:,0] = PI*(X[:,0]-0.5)/5        #latitude
X[:,1] = TWO_PI*(X[:,1]-0.5)/10    #longitude

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plot_sphere(ax)
x, y, z = to_cartesian(X[:,0], X[:,1])
ax.scatter3D(x, y, z)
ax.scatter3D([cos(pt_out[0])*cos(pt_out[1])], [sin(pt_out[0])], [cos(pt_out[0])*sin(pt_out[1])], c='g')
plt.show()

F = border(X, pt_out, min_dist, max_dist, True)

