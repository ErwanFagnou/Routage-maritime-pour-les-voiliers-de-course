from math import *
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import griddata                  #Ajustement du format des polaires (interpolation)
from skimage.transform import rescale, resize           #Ajustement du format des cartes de vent (agrandissement par un facteur)

from matplotlib.image import imread                     #Lecture de l'image contenant les obstacles
from scipy.io import netcdf                             #Lecture du format netcdf pour les courants
import json                                             #Lecture du format json pour les vents
import csv                                              #Lecture du format csv pour les polaires

from mpl_toolkits.mplot3d import Axes3D                 #Affichage 3D des points
import time


# Bugs dus surement aux erreurs de calcul (car grande influence de epsilon)
#   -> diminuer le nbre de points, donc augmenter leur distance / min_dist permet d'eviter les pts trop proches / croiser les doigts



# Vents: intervalle de 3h pendant 18 jours
# Courants: environ 24h (varie legerement) pendant 17 jours





# Données     ======================================================================
Rt = 3443.92                                            #Rayon de la Terre (milles marins)
PI = np.pi
TWO_PI = 2*PI
epsilon = 10**(-15)                                     #Marge d'erreur des calculs
                                                        #Saint-Malo: 48°38' N   2°02' W   ->   lat = 48.13°   long = -2.03°
d_lat = radians(48.63)%TWO_PI
d_long = radians(-2.03)%TWO_PI
                                                        #Guadeloupe: 16°13' N   61°32' W  ->   lat = 16.22°   long = -61.53°
a_lat = radians(16.22)%TWO_PI
a_long = radians(-61.53)%TWO_PI

m = 2                                                   #Facteur d'agrandissement de la grille
#m >= 2 !!!

start_time = 1541335500                                 #Heure de départ
dt = 10800                                              #Pas entre les itérations

n_points = 20                                           #Nombre max de points générés par un point à chaque itération
n_angle = int(360/n_points)                             #Pas entre les angles
min_dist = 1                                            #Distance minimale entre deux points
max_dist = 150                                          #Distance maximale entre deux points



# Polaires    ======================================================================
def load_polar(file):
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
def load_winds(file):
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

    currents_index = 0
    for k in range(0, len(data)):
        for i in range(0, len(np.array(data[k]['Datas']))):
            time = data[k]['Datas'][i]['timestamp']
            wind_time.append(time)
            while currents_index+1 < len(currents_time) and time > currents_time[currents_index+1]:
                currents_index += 1

            
            wind = np.array(data[k]['Datas'][i]['datas'], dtype=float)
            wind = np.reshape(wind, (97, 173, 2))
            wind[:, :, 1] = (450-wind[:, :, 1]*360/255)%360                                                     #Transformation angles pour sens trigo
            wind[:, :, 0] = wind[:, :, 0]/2.62

            u = rescale(wind[:, :, 0] * np.cos(np.deg2rad(wind[:, :, 1])), m, multichannel=False)             #On agrandit les projections sur Ox et Oy
            v = rescale(wind[:, :, 0] * np.sin(np.deg2rad(wind[:, :, 1])), m, multichannel=False)
            
            u2 = u - currents_data[currents_index, :, :, 0]
            v2 = v - currents_data[currents_index, :, :, 1]
            
            angle = (np.arctan2(v2, u2))%TWO_PI
            module = np.sqrt(u2**2 + v2**2)
            
            wind_data.append(np.swapaxes(np.swapaxes([module, angle], 0, 2), 0, 1))
    return np.array(wind_data), np.array(wind_time), min_long, max_long, min_lat, max_lat, m*nb_long, m*nb_lat

# Obstacles   ======================================================================
def load_obstacles(file):
    obstacles = imread(file)[:, :, 0]
    obstacles = resize(obstacles, (nb_lat, nb_long))>=0.5
    return obstacles

# Courants   ======================================================================
def load_currents(file):
    nc = netcdf.netcdf_file(file, 'r', mmap=False)
    s = nc.variables['uo'][:].shape
    f = int(m/2)
    
    currents_time = nc.variables['time'][:]*3600-631155600                  # Conversion en timestamp
    currents_data = np.zeros((s[0], f*(s[2]+1), f*(s[3]+1), 2))
    for i in range(len(currents_time)):
        u = np.flip(nc.variables['uo'][i][0], 1)
        v = np.flip(nc.variables['vo'][i][0], 1)
        u[u > 10**19] = 0
        v[v > 10**19] = 0
        
        u2 = rescale(u, f, multichannel=False)
        v2 = rescale(v, f, multichannel=False)
        
        currents_data[i, f:, f:, 0] = u2 * 1.9438612860586                  #Conversion en milles marin / h (noeuds)
        currents_data[i, f:, f:, 1] = v2 * 1.9438612860586

        currents_data[i, 0, f:, :] = currents_data[i, f, f:, :]             #On augmente la largeur et la longueur de 1 pour avoir les memes dimensions que wind_data
        currents_data[i, f:, 0, :] = currents_data[i, f:, f, :]
        currents_data[i, 0, 0, :] = currents_data[i, f, f, :]
    
    return currents_data, currents_time

# Point suivant     ======================================================================
def move_forward(lat, long, angle, time):
    y = int(nb_lat - 1 - (lat-min_lat)/((max_lat-min_lat)%TWO_PI)*(nb_lat-1))
    x = int((long-min_long)/((max_long-min_long)%TWO_PI)*(nb_long-1))
    
    w_index = int((time-wind_time[0])/(wind_time[1]-wind_time[0]))
    w_angle = wind_data[w_index, y, x, 1]
    w_speed = wind_data[w_index, y, x, 0]
    gap = abs(w_angle-angle)%TWO_PI
    
    v = polar[int(degrees(min(gap, TWO_PI-gap))),  int(w_speed)]
    
    s = dt * v / 3600 / Rt        #vitesse: milles/h   Rt: milles   -> donne la dist sur une sphere de rayon 1 (un angle donc)
    
    c_index = int((time-currents_time[0])/(currents_time[1]-currents_time[0]))
    c_lat = dt * currents_data[c_index, y, x, 0] / 3600 / Rt / cos(lat)
    c_long = dt * currents_data[c_index, y, x, 0] / 3600 / Rt

    
    lat2 = (asin(sin(lat)*cos(s) + cos(lat)*sin(s)*sin(angle))  +  c_lat)%TWO_PI
    long2 = (long  +  atan2(sin(s)*cos(angle), (cos(lat)*cos(s)-sin(lat)*sin(s)*sin(angle)))  +  c_long)%TWO_PI
    return lat2, long2

# Angle          ======================================================================
def angleGCR(lat1, long1, lat2, long2):
    long12 = long2-long1
##    p = sin(lat1+lat2) #tentative d'optimisation du calcul...
##    q = sin(lat1-lat2)
##    return (atan2((p-q)-(p+q)*cos(long12),2*cos(lat2)*sin(long12)))%TWO_PI
    c2 = cos(lat2)
    return (atan2(cos(lat1)*sin(lat2)-sin(lat1)*c2*cos(long12),c2*sin(long12)))%TWO_PI

def angleGCR2(A, B):
    return angleGCR(A[0], A[1], B[0], B[1])

# Distance    ======================================================================
def dist(lat1, long1, lat2, long2):
    long12 = long2-long1
    p = cos(lat1+lat2)
    q = cos(lat1-lat2)
    s = ((q-p)+(q+p)*cos(long12))/2
    if s>1:                 #si erreurs de calculs (ex: s=1.0000000000000002)
        s=1
    elif s<-1:
        s=-1
    return Rt*acos(s)

def dist2(A, B):
    return dist(A[0], A[1], B[0], B[1])


# Test sur intervalle circulaire =======================================================
def is_between_angles(a, a_min, a_max):
    a = (a-a_min)%TWO_PI
    a_max = (a_max-a_min)%TWO_PI
    return a <= a_max

# Coord spheriques -> cartesien ================================================================
def to_cartesian(lat, long):
    c = cos(lat)
    return c*cos(long), sin(lat), c*sin(long)

#Fonctions de vecteurs ========================================================================
def equal(A, B):
    return A[0] == B[0] and A[1] == B[1]

def cross(A, B):
    return [A[i-2]*B[i-1] - A[i-1]*B[i-2] for i in range(3)]

def dot(A, B):
    return A[0]*B[0]+A[1]*B[1]+A[2]*B[2]

def norme2(V):
    return V[0]**2 + V[1]**2 + V[2]**2

# Intersection d'arcs    ======================================================================
def do_arcs_intersect(A, B, C, D):  
    if equal(A, B) or equal(A, C) or equal(A, D) or equal(B, C) or equal(B, D) or equal(C, D):
        return False
    Ac = to_cartesian(A[0], A[1])
    Bc = to_cartesian(B[0], B[1])
    Cc = to_cartesian(C[0], C[1])
    Dc = to_cartesian(D[0], D[1])
    
    V2 = cross(Cc, Dc)
    if dot(V2, Ac)*dot(V2, Bc) < -epsilon: #ie A et B de chaque coté du plan OCD (si =0 on ne compte pas l'intersection)
        V1 = cross(Ac, Bc)
        if dot(V1, Cc)*dot(V1, Dc) < -epsilon and norme2(cross(V1, V2)) > epsilon: #idem pour C et D, et les points ne sont pas alignés
            V3 = cross(Ac, Cc)
            if dot(V3, Bc)*dot(V3, Dc) > epsilon: #ie si B et D appartiennent au meme hemisphere (A et C dedans car le définissent)
                return True
    return False

# Calcul de la bordure du nuage de pts ======================================================================
def get_border(X, pt_out, min_dist, max_dist):
    best_dist = dist2(X[0], pt_out)
    first_pt = X[0]
    for P in X:
        d = dist2(P, pt_out)
        if d < best_dist:
            best_dist = d
            first_pt = P
    F = [first_pt]
    X_mask =  np.ones((len(X),))
    prev_angle = angleGCR2(first_pt, pt_out) % TWO_PI
    while len(F)==1 or ((not equal(F[-1],F[0])) and len(F)==2) or ((not ( equal(F[-1], F[1]) and equal(F[-2], F[0]) ) ) and len(F)>=3):

        if len(F)>1000:
            break
        
        best_angle = TWO_PI + 1
        best_dist = max_dist+1
        next_pt = F[-1]
        for p in range(len(X)):
            if X_mask[p]==1:
                P = X[p]
                d = dist2(F[-1], P)
                if (not equal(P, F[-1])) and d <= max_dist:
                    angle = (angleGCR2(F[-1], P) - prev_angle) % TWO_PI
                    if angle <= 0:
                        angle = TWO_PI
                    if min_dist <= d:
                        if angle < best_angle or (angle==best_angle and d < best_dist):
                            intersect = do_arcs_intersect(F[-1], P, pt_out, first_pt)
                            i = 0
                            while not intersect and i < len(F)-1:
                                intersect = do_arcs_intersect(F[-1], P, F[i], F[i+1])
                                i+=1
                            if not intersect:
                                best_angle = angle
                                best_dist = d
                                next_pt = P
                    else:
                        X_mask[p] = 0

        prev_angle = angleGCR2(next_pt, F[-1])
        F.append(next_pt)
    if len(F)>1000:
        print("Erreur!!!")
        X2 = np.array(X)[X_mask==1]
        F2 = np.array(F)
        plt.scatter(X2[:,1], X2[:,0])
        plt.scatter(F2[:,1], F2[:,0], c='r')

        plt.plot(F2[:,1], F2[:,0], 'r')
        plt.show()
    if len(F)<=3:
        return F[:-1]
    else:
        return F[:-2]








# Algorithme   =====================================================================

# Polaires
polar = np.zeros((181, 71))
polarFiles = ['vpp_1_1.csv', 'vpp_1_2.csv', 'vpp_1_4.csv', 'vpp_1_8.csv', 'vpp_1_16.csv', 'vpp_1_32.csv', 'vpp_1_64.csv']
for file in polarFiles:
    polar = np.maximum(polar, load_polar('data/polaires/'+file))
print('- Polaires chargées')
##plt.imshow(np.repeat(np.expand_dims(polar/np.max(polar), 2), 3, axis=2))
##plt.show()


# Courants
currents_data, currents_time = load_currents('data/currents.nc')
print('- Carte des courants chargée')
##for i in range(len(currents_time)):
##    print(np.max(currents_data[i, :, :, 0]))
##    plt.imshow(np.repeat(np.expand_dims((currents_data[i, :, :, 0]-np.min(currents_data[i, :, :, 0]))/(np.max(currents_data[i, :, :, 0]) - np.min(currents_data[i, :, :, 0])), 2), 3, axis=2))
##    plt.show()
##    plt.imshow(np.repeat(np.expand_dims((currents_data[i, :, :, 1]-np.min(currents_data[i, :, :, 1]))/(np.max(currents_data[i, :, :, 1]) - np.min(currents_data[i, :, :, 1])), 2), 3, axis=2))
##    plt.show()

# Vent
wind_data, wind_time, min_long, max_long, min_lat, max_lat, nb_long, nb_lat = load_winds("data/wind.json")
print('- Carte des vents chargée')
##print(min_lat, min_long, max_lat, max_long, nb_lat, nb_long)
##
##for i in range(len(wind_time)):
##    print(np.max(wind_data[i, :, :, 0]))
##    plt.imshow(np.repeat(np.expand_dims(wind_data[i, :, :, 0]/np.max(wind_data[i, :, :, 0]), 2), 3, axis=2))
##    plt.show()
##    plt.imshow(np.repeat(np.expand_dims(wind_data[i, :, :, 1]/np.max(wind_data[i, :, :, 1]), 2), 3, axis=2))
##    plt.show()



# Obstacles
obstacles = load_obstacles('data/obstacles_north_atlantic.png')
print('- Carte des obstacles chargée')
##plt.imshow(np.repeat(np.expand_dims(obstacles/np.max(obstacles), 2), 3, axis=2))
##plt.show()




depart = [d_lat, d_long]
iso = [[[d_lat, d_long, 0]]]
t0 = time.time()
for i in range(79):
    print(i, "Temps:", i*dt, end=" ")
    new_iso = []
    n = len(iso[i])
    for k in range(len(iso[i])):
        lat = iso[i][k][0]
        long = iso[i][k][1]
        blocked = False
        if obstacles[int(nb_lat - 1 - (lat-min_lat)/((max_lat-min_lat)%TWO_PI)*(nb_lat-1)) , int((long-min_long)/((max_long-min_long)%TWO_PI)*(nb_long-1))]:
            previous_angle = int(degrees(angleGCR2(iso[i][k], iso[i][k-1])))
            next_angle = int(degrees(angleGCR2(iso[i][k], iso[i][(k+1)%n])))
            if next_angle <= previous_angle:
                next_angle += 360
            for angle in range(previous_angle, next_angle + n_angle, n_angle):
                lat2, long2 = move_forward(lat, long, radians(angle), i*dt+start_time)
                y2 = int(nb_lat - 1 - (lat2-min_lat)/((max_lat-min_lat)%TWO_PI)*(nb_lat-1))
                x2 = int((long2-min_long)/((max_long-min_long)%TWO_PI)*(nb_long-1))
                if is_between_angles(lat2, min_lat, max_lat) and is_between_angles(long2, min_long, max_long) and obstacles[y2, x2]:
                    new_iso.append([lat2, long2, k])
                else:
                    blocked = True
        else:
            blocked = True
        if blocked:
            new_iso.append([lat, long, k])
    #new_iso = np.unique(new_iso, axis=0).tolist()
    B = get_border(new_iso, [a_lat, a_long], min_dist, max_dist)
    print("n:", len(new_iso), " h:", len(B))
    iso.append(B)



    
print("Temps:", time.time()-t0)



plt.imshow(np.repeat(np.expand_dims(obstacles/np.max(obstacles), 2), 3, axis=2))
for lst in iso:
    X = []
    Y = []
    lst = np.append(lst, [lst[0]], axis=0)
    for point in lst:
        y = nb_lat - 1 - (point[0]-min_lat)/((max_lat-min_lat)%TWO_PI)*(nb_lat-1)-0.5
        x = (point[1]-min_long)/((max_long-min_long)%TWO_PI)*(nb_long-1)-0.5
        X.append(x)
        Y.append(y)
    plt.plot(X, Y, 'r')
plt.title("Isochrones (dt=" + str(dt) + ", n_pt=" + str(n_points) + ", min_dist=" + str(min_dist) + ", max_dist=" + str(max_dist) + ")")
plt.show()
