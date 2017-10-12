import csv
import math
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

wthetas = np.zeros(360)
#normal vector
nv = np.array([0,1,0],dtype=float)
cp = np.array([1,0,0],dtype=float)

with open('final_vessels.csv') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        #print(row['x1'])
        p0 = [row['x0'], row['y0'], row['z0']]
        p1 = [row['x1'], row['y1'], row['z1']]
        p0=np.array(p0,dtype=float)
        p1=np.array(p1,dtype=float)
        length = np.linalg.norm(p1-p0)
        assert(length >= 0.0)
        dp0 = np.dot(nv, p1 -p0)/length
        dp1 = np.dot(cp, p1-p0)/length
        
        theta = np.degrees(np.arccos(dp0))
        ftheta = np.degrees(np.arccos(dp1)) 
        if theta > 90:
            ftheta = 360 - ftheta
        ftheta = math.floor(ftheta)
        wthetas[ftheta] += length


y_pos = np.arange(len(wthetas))
plt.bar(y_pos, wthetas, align='center', alpha=0.5)
#plt.xticks(y_pos, objects)
plt.ylabel('Length')
plt.title('Angle')
 
plt.show()