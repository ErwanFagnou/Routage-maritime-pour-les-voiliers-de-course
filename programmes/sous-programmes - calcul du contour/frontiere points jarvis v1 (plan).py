import random
import matplotlib.pyplot as plt
import numpy as np
from math import *
import time

PI = 3.14159265
TWO_PI = 6.28318531

def dist(A, B):
    return sqrt((A[0]-B[0])**2 + (A[1]-B[1])**2)

def angle_rot(A, B):
    return atan2(B[1]-A[1], B[0]-A[0])%TWO_PI

def is_between_angles(a, a_min, a_max):
    if a_min <= a_max:
        return a_min < a and a < a_max
    else:
        return a_min < a or a < a_max

def do_segments_intersect(A, B, C, D):
    a1 = angle_rot(A, B)
    a2 = angle_rot(C, D)
    if (a1-a2)%TWO_PI <= PI:
        return is_between_angles(a1, angle_rot(A, D), angle_rot(A, C)) and is_between_angles(a2, angle_rot(C, A), angle_rot(C, B)) and is_between_angles(angle_rot(B, A) , angle_rot(B, C), angle_rot(B, D)) and is_between_angles(angle_rot(D, C), angle_rot(D, B), angle_rot(D, A))
    else:
        return is_between_angles(a1, angle_rot(A, C), angle_rot(A, D)) and is_between_angles(a2, angle_rot(C, B), angle_rot(C, A)) and is_between_angles(angle_rot(B, A) , angle_rot(B, D), angle_rot(B, C)) and is_between_angles(angle_rot(D, C), angle_rot(D, A), angle_rot(D, B))

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
    if show:
        F2 = np.array(F)
        plt.scatter(X[:,0], X[:,1])#, marker='x')
        plt.scatter([pt_out[0]], [pt_out[1]], c='r')
        plt.plot([pt_out[0], first_pt[0]], [pt_out[1], first_pt[1]], 'g')
        plt.plot(F2[:,0], F2[:,1], 'r')
        plt.show()
    return F




n = 1000
max_dist = 1.5
pt_out = [5,5]


X = np.random.normal(0,1,(n,2))
#X = (2*X-1)/((0.5+np.abs(2*X-1)**3))
plt.scatter(X[:,0], X[:,1])#, marker='x')
plt.scatter([pt_out[0]], [pt_out[1]], c='r')
plt.show()

F = border(X, pt_out, 0, max_dist, show=True)

