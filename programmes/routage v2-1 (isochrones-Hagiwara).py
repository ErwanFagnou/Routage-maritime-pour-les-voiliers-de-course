from math import *
import json
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata                  #Ajustement du format des polaires (interpolation)
from skimage.transform import rescale, resize           #Ajustement du format des cartes de vent (agrandissement par un facteur)
from matplotlib.image import imread                     #Lecture de l'image contenant les obstacles

from mpl_toolkits.mplot3d import Axes3D                 #Affichage 3D des points

import time as tm

# Données     ======================================================================
Rt = 3443.92 #miles marins
PI = np.pi
TWO_PI = 2*PI

    #Saint-Malo: 48°38' N   2°02' W   ->   lat = 48.13   long = -2.03
d_lat = radians(48.63)%TWO_PI
d_long = radians(-2.03)%TWO_PI
    #Guadeloupe: 16°13' N   61°32' W  ->   lat = 16.22   long = -61.53
a_lat = radians(16.22)%TWO_PI
a_long = radians(-61.53)%TWO_PI

m = 4                                                                                       #Facteur d'agrandissement de la grille

start_time = 1541335500
dt = 10800

n_segments = 500                                                                          #Nombre de segments découpant les isochrones
n_points = 30



# ==================================================================================















# Polaires    ======================================================================
def loadPolar(file):
    reader = csv.reader(open(file, newline=''), delimiter=';')
    points = []
    values = []
    for row in reader:
        if row[0] == 'TWA':
            first = [int(i) for i in row[1:]]
        else:
            for i in range(len(row)-1):
                points.append([int(row[0]), first[i]])
                values.append(float(row[i+1]))
    grid_x, grid_y = np.mgrid[0:181, 0:71]
    return griddata(points, values, (grid_x, grid_y), method='cubic')

# Vents       ======================================================================
def loadWinds(file):
    f = open(file, "r")
    data = json.loads(f.read())
    wind_data = []
    wind_time = []
    min_long = radians(data[0]['Datas'][0]['minlong'])%TWO_PI
    max_long = radians(data[0]['Datas'][0]['maxlong'])%TWO_PI
    min_lat = radians(data[0]['Datas'][0]['minlat'])%TWO_PI
    max_lat = radians(data[0]['Datas'][0]['maxlat'])%TWO_PI
    nb_long = data[0]['Datas'][0]['nblong']
    nb_lat = data[0]['Datas'][0]['nblat']
    for k in range(0, len(data)):
        for i in range(0, len(np.array(data[k]['Datas']))):
            wind_time.append(np.array(data[k]['Datas'][i]['timestamp']))
            wind = np.array(data[k]['Datas'][i]['datas'], dtype=float)
            wind = np.reshape(wind, (97, 173, 2))
            wind[:, :, 1] = (450-wind[:, :, 1]*360/255)%360                                                                             #Transformation angles pour sens trigo
            
            cosAngle = rescale(np.cos(np.deg2rad(wind[:, :, 1]))/2+0.5, m, multichannel=False)*2-1
            sinAngle = rescale(np.sin(np.deg2rad(wind[:, :, 1]))/2+0.5, m, multichannel=False)*2-1
            angle = (np.arctan2(sinAngle, cosAngle))%TWO_PI
            
            speed = rescale(wind[:, :, 0]/255, m, multichannel=False)*255/2.62
            wind_data.append(np.swapaxes(np.swapaxes([speed, angle], 0, 2), 0, 1))
    return np.array(wind_data), np.array(wind_time), min_long, max_long, min_lat, max_lat, m*nb_long, m*nb_lat

# Obstacles   ======================================================================
def loadObstacles(file):
    obstacles = imread(file)[:, :, 0]
    obstacles = resize(obstacles, (nb_lat, nb_long))>=0.5
    return obstacles

# Vitesse     ======================================================================
def speed(lat, long, time, angle):
    #print(lat, long, min_lat, min_long, angle)
    y = int(nb_lat - 1 - (lat-min_lat)/((max_lat-min_lat)%TWO_PI)*(nb_lat-1))
    x = int((long-min_long)/((max_long-min_long)%TWO_PI)*(nb_long-1))
    #print(x, y)

    if obstacles[y, x]:
        w_index = int((time-wind_time[0])/(wind_time[1]-wind_time[0]))
        w_angle = wind_data[w_index, y, x, 1]
        w_speed = wind_data[w_index, y, x, 0]
        gap = abs(w_angle-angle)%TWO_PI
        return polar[int(degrees(min(gap, TWO_PI-gap))),  int(w_speed)]
    else:
        return 0

# Point suivant     ======================================================================
def moveForward(lat, long, angle, v):
    s = dt*v/3600/Rt        #vitesse: miles/h   Rt: miles   -> donne la dist sur une sphere de rayon 1
    lat2 = asin(  sin(lat)*cos(s) + cos(lat)*sin(s)*sin(angle) )%TWO_PI
    long2 = (long + atan2(sin(s)*cos(angle), (cos(lat)*cos(s)-sin(lat)*sin(s)*sin(angle))))%TWO_PI
    return lat2, long2

# Point suivant     ======================================================================
def angleGCR(lat1, long1, lat2, long2):
    long12 = long2-long1
    
    return (PI/2-atan2(cos(lat2)*sin(long12),(cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(long12))))%TWO_PI
    
# Distance    ======================================================================
def dist(lat1, long1, lat2, long2):
    long12 = long2-long1
    
    return 60*180/PI*acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(long12)) # 1 degré = 60 miles marins = pi/180 radians
# ==================================================================================

def is_between_angles(a, a_min, a_max):
    if a_min <= a_max:
        return a_min < a and a < a_max
    else:
        return a_min < a or a < a_max


















# Algorithme   =====================================================================

# Polaires
polar = np.zeros((181, 71))
polarFiles = ['vpp_1_1.csv', 'vpp_1_2.csv', 'vpp_1_4.csv', 'vpp_1_8.csv', 'vpp_1_16.csv', 'vpp_1_32.csv', 'vpp_1_64.csv']
for file in polarFiles:
    polar = np.maximum(polar, loadPolar('data/polaires/'+file))
print('- Polaires chargées')
#plt.imshow(np.repeat(np.expand_dims(polar/np.max(polar), 2), 3, axis=2))
#plt.show()

# Vent
wind_data, wind_time, min_long, max_long, min_lat, max_lat, nb_long, nb_lat = loadWinds("data/wind.json")
print('- Carte des vents chargée')

# Obstacles
obstacles = loadObstacles('data/obstacles_north_atlantic.png')
print('- Carte des obstacles chargée')
#plt.imshow(np.repeat(np.expand_dims(obstacles/np.max(obstacles), 2), 3, axis=2))
#plt.show()


t0 = tm.time()


depart = [d_lat, d_long]

iso = [[[d_lat, d_long, 0]]]
parents = [[]]
for i in range(85):
    print(i, "Time:", i*dt)
    segments = [[] for n in range(n_segments)]
    for k in range(len(iso[i])):
        lat = iso[i][k][0]
        long = iso[i][k][1]
        for angle in range(0, 360, int(360/n_points)):
            v = speed(lat, long, i*dt+start_time, radians(angle))
            lat2, long2 = moveForward(lat, long, radians(angle), v)
            if is_between_angles(lat2, min_lat, max_lat) and is_between_angles(long2, min_long, max_long):
                #print("new point: ", lat2, long2)
                angleDepart = angleGCR(d_lat, d_long, lat2, long2)
                segments[int(angleDepart*n_segments/TWO_PI)].append([lat2, long2, k])
                #print(lat2, long2, angleDepart*n_segments/360)
    new_iso = []
    for s in segments:
        if s!=[]:
            best = dist(d_lat, d_long, s[0][0], s[0][1])
            index = 0
            for i in range(1, len(s)):
                d = dist(d_lat, d_long, s[i][0], s[i][1])
                if d > best:
                    best = d
                    index = i
            new_iso.append(s[index])
    iso.append(new_iso)
##    
##    plt.imshow(obstacles)
##    X = []
##    Y = []
##    new_iso.append(new_iso[0])
##    for point in new_iso:
##        y = nb_lat - (point[0]-min_lat)/(max_lat-min_lat)*nb_lat
##        x = (point[1]-min_long)/(max_long-min_long)*nb_long
##        X.append(x)
##        Y.append(y)
##    plt.plot(X, Y, 'r')
##    plt.show()


print(tm.time()-t0)

plt.imshow(np.repeat(np.expand_dims(obstacles/np.max(obstacles), 2), 3, axis=2))
for lst in iso:
    X = []
    Y = []
    lst.append(lst[0])
    for point in lst:
        y = nb_lat - (point[0]-min_lat)/((max_lat-min_lat)%TWO_PI)*(nb_lat-1)
        x = (point[1]-min_long)/((max_long-min_long)%TWO_PI)*(nb_long-1)
        X.append(x)
        Y.append(y)
    plt.plot(X, Y, 'r')
plt.title("Isochrones (dt=" + str(dt) + ", n_seg=" + str(n_segments) + ", n_pt=" + str(n_points))
plt.show()



dt = 10800

n_segments = 500                                                                          #Nombre de segments découpant les isochrones
n_points = 50


























