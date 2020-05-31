import json
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata                  #Ajustement du format des polaires (interpolation)
from math import *
from skimage.transform import rescale, resize           #Ajustement du format des cartes de vent (agrandissement par un facteur)
from matplotlib.image import imread                     #Lecture de l'image contenant les obstacles
import time as tm

"""
Le programme peut être un peu difficile à lire car j'ai fait mon maximum pour
l'optimiser et réduire le temps de calcul. Par exemple j'ai imposé une vitesse
minimale, et je ne m'intéresse pas à une case que j'ai déjà traité il y a
plus de forget_steps tours.
"""

# Polaires    =====================================================================
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
    min_long = data[0]['Datas'][0]['minlong']
    max_long = data[0]['Datas'][0]['maxlong']
    min_lat = data[0]['Datas'][0]['minlat']
    max_lat = data[0]['Datas'][0]['maxlat']
    nb_long = data[0]['Datas'][0]['nblong']
    nb_lat = data[0]['Datas'][0]['nblat']
    for k in range(0, len(data)):
        for i in range(0, len(np.array(data[k]['Datas']))):
            wind_time.append(np.array(data[k]['Datas'][i]['timestamp']))
            wind = np.array(data[k]['Datas'][i]['datas'], dtype=float)
            wind = np.reshape(wind, (97, 173, 2))
            wind[:, :, 1] = (450-wind[:, :, 1]*360/255)%360         #Transformation angles pour sens trigo
            
            cosAngle = rescale(np.cos(np.deg2rad(wind[:, :, 1]))/2+0.5, m, multichannel=False)*2-1
            sinAngle = rescale(np.sin(np.deg2rad(wind[:, :, 1]))/2+0.5, m, multichannel=False)*2-1
            angle = (np.rad2deg(np.arctan2(sinAngle, cosAngle))+360)%360
            
            speed = rescale(wind[:, :, 0]/255, m, multichannel=False)*255/2.61
            wind_data.append(np.swapaxes(np.swapaxes([speed, angle], 0, 2), 0, 1))
    return np.array(wind_data), np.array(wind_time), min_long, max_long, min_lat, max_lat, m*nb_long, m*nb_lat

# Obstacles   ======================================================================
def loadObstacles(file):
    obstacles = imread(file)[:, :, 0]
    obstacles = resize(obstacles, (nb_lat, nb_long))>=0.5
    #print(obstacles.shape)
    return obstacles

# Données     ======================================================================

    #Saint-Malo: 48°38' N   2°02' W   ->   lat = 48.13   long = -2.03
d_lat = 48.13
d_long = -2.03
    #Guadeloupe: 16°13' N   61°32' W  ->   lat = 16.22   long = -61.53
a_lat = 16.22
a_long = -61.53

m = 1                               #Facteur d'agrandissement de la grille (gros impact sur les performances, dépendance en m^2)

v_min = 0.1                         #Vitesse minimale à prendre en compte (évite les divisions par zéro)
forget_steps = 5                    #Nombre d'itérations avant de ne plus utiliser une cellule de la grille


# Algorithme   =====================================================================

# Polaires
polar = np.zeros((181, 71))
polarFiles = ['vpp_1_1.csv', 'vpp_1_2.csv', 'vpp_1_4.csv', 'vpp_1_8.csv', 'vpp_1_16.csv', 'vpp_1_32.csv', 'vpp_1_64.csv']
for file in polarFiles:
    polar = np.maximum(polar, loadPolar('data/polaires/'+file))
print('- Polaires chargées')
##plt.imshow(np.repeat(np.expand_dims(polar/np.max(polar), 2), 3, axis=2))
##plt.show()

# Vent
wind_data, wind_time, min_long, max_long, min_lat, max_lat, nb_long, nb_lat = loadWinds("data/wind.json")
print('- Carte des vents chargée')

# Obstacles
obstacles = loadObstacles('data/obstacles_north_atlantic.png')
print('- Carte des obstacles chargée')
##plt.imshow(np.repeat(np.expand_dims(obstacles/np.max(obstacles), 2), 3, axis=2))
##plt.show()


t0 = tm.time()

depart = [nb_lat - int((d_lat-min_lat)/(max_lat-min_lat)*nb_lat), int((d_long-min_long)/(max_long-min_long)*nb_long)]
arrivee = [nb_lat - int((a_lat-min_lat)/(max_lat-min_lat)*nb_lat), int((a_long-min_long)/(max_long-min_long)*nb_long)]

state = np.zeros((nb_lat, nb_long), dtype=int)
parents = -np.ones((nb_lat, nb_long, 2), dtype=int)
time = -np.ones((nb_lat, nb_long), dtype=int)

#0: non visité,    1: déjà visité,    2: déjà visté et exploité
state[depart[0], depart[1]] = 1
time[depart[0], depart[1]] = 1541335500
timestep = 10800
actual_time = time[depart[0], depart[1]]
wind_index = 0

for n in range(1000):
    actual_time += timestep
    while wind_time[wind_index]+timestep < actual_time:
        wind_index += 1
    print(n,'Time:', actual_time-time[depart[0], depart[1]])
    while 1 in state:
        for i in range(0, nb_lat):
            for j in range(0, nb_long):
                if state[i, j] == 1:
                    state[i, j] = 2
                    
                    w_angle = int(wind_data[wind_index, i, j, 1])
                    w_speed = wind_data[wind_index, i, j, 0]
                    lat =  min_lat + i*(max_lat-min_lat)/nb_lat
#Sud
                    gap = abs(w_angle-270)
                    v = polar[min(gap, 360-gap),  int(w_speed)]
                    if v>v_min:
                        new_time = time[i, j] + 34.6/m/v*3600
                        if new_time<actual_time and i<nb_lat-1 and obstacles[i+1, j] and (state[i+1, j] == 0 or time[i+1, j] > new_time):
                            state[i+1, j] = 1+(state[i+1, j] == 2)
                            parents[i+1, j] = [i, j]
                            time[i+1, j] = new_time

#Nord
                    gap = abs(w_angle-90)
                    v = polar[min(gap, 360-gap),  int(w_speed)]
                    if v>v_min:
                        new_time = time[i, j] + 34.6/m/v*3600
                        if new_time<actual_time and i>0 and obstacles[i-1, j] and (state[i-1, j] == 0 or time[i-1, j] > new_time):
                            state[i-1, j] = 1+(state[i-1, j] == 2)
                            parents[i-1, j] = [i, j]
                            time[i-1, j] = new_time
#Est
                    gap = abs(w_angle)
                    v = polar[min(gap, 360-gap),  int(w_speed)]
                    if v>v_min:
                        new_time = time[i, j] + abs(cos(radians(lat)))*34.6/m/v*3600
                        if new_time<actual_time and j<nb_long-1 and obstacles[i, j+1] and (state[i, j+1] == 0 or time[i, j+1] > new_time):
                            state[i, j+1] = 1+(state[i, j+1] == 2)
                            parents[i, j+1] = [i, j]
                            time[i, j+1] = new_time
#Ouest
                    gap = abs(w_angle-180)
                    v = polar[min(gap, 360-gap),  int(w_speed)]
                    if v>v_min:
                        new_time = time[i, j] + abs(cos(radians(lat)))*34.6/m/v*3600
                        if new_time<actual_time and j>0 and obstacles[i, j-1] and (state[i, j-1] == 0 or time[i, j-1] > new_time):
                            state[i, j-1] = 1+(state[i, j-1] == 2)
                            parents[i, j-1] = [i, j]
                            time[i, j-1] = new_time
#Sud-est
                    gap = abs(w_angle-315)
                    v = polar[min(gap, 360-gap),  int(w_speed)]
                    if v>v_min:
                        new_time = time[i, j] + 3959*acos(sin(radians(lat))*sin(radians(lat-0.5/m))+cos(radians(lat))*cos(radians(lat-0.5/m))*cos(radians(0.5/m)))/v*3600
                        if new_time<actual_time and i<nb_lat-1 and j<nb_long-1 and obstacles[i+1, j+1] and (state[i+1, j+1] == 0 or time[i+1, j+1] > new_time):
                            state[i+1, j+1] = 1+(state[i+1, j+1] == 2)
                            parents[i+1, j+1] = [i, j]
                            time[i+1, j+1] = new_time
#Nord-est
                    gap = abs(w_angle-45)
                    v = polar[min(gap, 360-gap),  int(w_speed)]
                    if v>v_min:
                        new_time = time[i, j] + 3959*acos(sin(radians(lat))*sin(radians(lat+0.5/m))+cos(radians(lat))*cos(radians(lat+0.5/m))*cos(radians(0.5/m)))/v*3600
                        if new_time<actual_time and i>0 and j<nb_long-1 and obstacles[i-1, j+1] and (state[i-1, j+1] == 0 or time[i-1, j+1] > new_time):
                            state[i-1, j+1] = 1+(state[i-1, j+1] == 2)
                            parents[i-1, j+1] = [i, j]
                            time[i-1, j+1] = new_time
#Sud-ouest
                    gap = abs(w_angle-225)
                    v = polar[min(gap, 360-gap),  int(w_speed)]
                    if v>v_min:
                        new_time = time[i, j] + 3959*acos(sin(radians(lat))*sin(radians(lat-0.5/m))+cos(radians(lat))*cos(radians(lat-0.5/m))*cos(radians(0.5/m)))/v*3600
                        if new_time<actual_time and i<nb_lat-1 and j>0 and obstacles[i+1, j-1] and (state[i+1, j-1] == 0 or time[i+1, j-1] > new_time):
                            state[i+1, j-1] = 1+(state[i+1, j-1] == 2)
                            parents[i+1, j-1] = [i, j]
                            time[i+1, j-1] = new_time
#Nord-ouest
                    gap = abs(w_angle-135)
                    v = polar[min(gap, 360-gap),  int(w_speed)]
                    if v>v_min:
                        new_time = time[i, j] + 3959*acos(sin(radians(lat))*sin(radians(lat+0.5/m))+cos(radians(lat))*cos(radians(lat+0.5/m))*cos(radians(0.5/m)))/v*3600
                        if new_time<actual_time and i>0 and j>0 and obstacles[i-1, j-1] and (state[i-1, j-1] == 0 or time[i-1, j-1] > new_time):
                            state[i-1, j-1] = 1+(state[i-1, j-1] == 2)
                            parents[i-1, j-1] = [i, j]
                            time[i-1, j-1] = new_time


#Nord-ouest-ouest
                    gap = int(abs(w_angle-157.5))
                    v = polar[min(gap, 360-gap),  int(w_speed)]
                    if v>v_min:
                        new_time = time[i, j] + 3959*acos(sin(radians(lat))*sin(radians(lat+0.5/m))+cos(radians(lat))*cos(radians(lat+0.5/m))*cos(radians(1/m)))/v*3600
                        if new_time<actual_time and i>0 and j>1 and obstacles[i-1, j-2] and (state[i-1, j-2] == 0 or time[i-1, j-2] > new_time):
                            state[i-1, j-2] = 1+(state[i-1, j-2] == 2)
                            parents[i-1, j-2] = [i, j]
                            time[i-1, j-2] = new_time


#Sud-ouest-ouest
                    gap = int(abs(w_angle-202.5))
                    v = polar[min(gap, 360-gap),  int(w_speed)]
                    if v>v_min:
                        new_time = time[i, j] + 3959*acos(sin(radians(lat))*sin(radians(lat-0.5/m))+cos(radians(lat))*cos(radians(lat-0.5/m))*cos(radians(1/m)))/v*3600
                        if new_time<actual_time and i<nb_lat-1 and j>1 and obstacles[i+1, j-2] and (state[i+1, j-2] == 0 or time[i+1, j-2] > new_time):
                            state[i+1, j-2] = 1+(state[i+1, j-2] == 2)
                            parents[i+1, j-2] = [i, j]
                            time[i+1, j-2] = new_time

#Sud-sud-ouest
                    gap = int(abs(w_angle-247.5))
                    v = polar[min(gap, 360-gap),  int(w_speed)]
                    if v>v_min:
                        new_time = time[i, j] + 3959*acos(sin(radians(lat))*sin(radians(lat-1/m))+cos(radians(lat))*cos(radians(lat-1/m))*cos(radians(0.5/m)))/v*3600
                        if new_time<actual_time and i<nb_lat-2 and j>0 and obstacles[i+2, j-1] and (state[i+2, j-1] == 0 or time[i+2, j-1] > new_time):
                            state[i+2, j-1] = 1+(state[i+2, j-1] == 2)
                            parents[i+2, j-1] = [i, j]
                            time[i+2, j-1] = new_time

#Sud-sud-est
                    gap = int(abs(w_angle-292.5))
                    v = polar[min(gap, 360-gap),  int(w_speed)]
                    if v>v_min:
                        new_time = time[i, j] + 3959*acos(sin(radians(lat))*sin(radians(lat-1/m))+cos(radians(lat))*cos(radians(lat-1/m))*cos(radians(0.5/m)))/v*3600
                        if new_time<actual_time and i<nb_lat-2 and j<nb_long-1 and obstacles[i+2, j+1] and (state[i+2, j+1] == 0 or time[i+2, j+1] > new_time):
                            state[i+2, j+1] = 1+(state[i+2, j+1] == 2)
                            parents[i+2, j+1] = [i, j]
                            time[i+2, j+1] = new_time
    
    state[(state==2) & (time>actual_time-forget_steps*timestep)] = 1
    if state[arrivee[0], arrivee[1]] == 1 or wind_index==151:
        break
    
if state[arrivee[0], arrivee[1]] == 0:
    print('Error: no path found')
    plt.imshow(state*100)
    plt.show()
else:
    print('Path found')
    trajet = np.array([arrivee], dtype=int)
    angle = []
    force = []
    while trajet[-1, 0] != depart[0] or trajet[-1, 1] != depart[1]:
        parent = parents[trajet[-1, 0], trajet[-1, 1]]
        trajet = np.append(trajet, [[parent[0], parent[1]]], axis=0)

    print(tm.time()-t0)

    trajetX = trajet[:, 0]
    trajetY = trajet[:, 1]

    plt.imshow(np.repeat(np.expand_dims(obstacles/np.max(obstacles), 2), 3, axis=2))
    plt.plot(trajetY, trajetX, 'r')
    plt.show()
