from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import healpy as hp

fig = plt.figure(figsize=(10,10))

ax = fig.add_subplot(111, projection='3d')
ax.set_aspect("equal")

cs_lon = 207.8
cs_lat = -56.3
cs_vec = hp.pixelfunc.ang2vec(cs_lon, cs_lat, lonlat = True)
cs_vec_norm = (1/np.sqrt(cs_vec[0]**2+cs_vec[1]**2+cs_vec[2]**2))*cs_vec

dom_1_lon = 117.8
dom_1_lat = 0
dom_1_vec = hp.pixelfunc.ang2vec(dom_1_lon, dom_1_lat, lonlat = True)
dom_1_vec_norm = (1/np.sqrt(dom_1_vec[0]**2+dom_1_vec[1]**2+dom_1_vec[2]**2))*dom_1_vec
dom_2_lon = 207.8
dom_2_lat = 33.7
dom_2_vec = hp.pixelfunc.ang2vec(dom_2_lon, dom_2_lat, lonlat = True)
dom_2_vec_norm = (1/np.sqrt(dom_2_vec[0]**2+dom_2_vec[1]**2+dom_2_vec[2]**2))*dom_2_vec

ax.plot(np.array((dom_1_vec_norm[0], -dom_1_vec_norm[0])),np.array((dom_1_vec_norm[1], -dom_1_vec_norm[1])),np.array((dom_1_vec_norm[2], -dom_1_vec_norm[2])), color = 'k')
ax.plot(np.array((dom_2_vec_norm[0], -dom_2_vec_norm[0])),np.array((dom_2_vec_norm[1], -dom_2_vec_norm[1])),np.array((dom_2_vec_norm[2], -dom_2_vec_norm[2])), color = 'k')
ax.plot(np.array((0,cs_vec_norm[0])), np.array((0,cs_vec_norm[1])), np.array((0,cs_vec_norm[2])), color = 'g')

r = 1
t = np.arange(0,2*np.pi+2*np.pi/360, 2*np.pi/360)
cos_t = np.cos(t)
sin_t = np.sin(t)

p_1_a = dom_2_vec_norm*np.array((0.1/np.sqrt(2),0.05, 0.05))
circle_1_a = np.zeros((3, len(t)))
circle_1_a[0] = p_1_a[0] + r*cos_t*dom_1_vec_norm[0]+r*sin_t*cs_vec_norm[0]
circle_1_a[1] = p_1_a[1] + r*cos_t*dom_1_vec_norm[1]+r*sin_t*cs_vec_norm[1]
circle_1_a[2] = p_1_a[2] + r*cos_t*dom_1_vec_norm[2]+r*sin_t*cs_vec_norm[2]

p_1_b = dom_2_vec_norm*-np.array((0.1/np.sqrt(2),0.05, 0.05))
circle_1_b = np.zeros((3, len(t)))
circle_1_b[0] = p_1_b[0] + r*cos_t*dom_1_vec_norm[0]+r*sin_t*cs_vec_norm[0]
circle_1_b[1] = p_1_b[1] + r*cos_t*dom_1_vec_norm[1]+r*sin_t*cs_vec_norm[1]
circle_1_b[2] = p_1_b[2] + r*cos_t*dom_1_vec_norm[2]+r*sin_t*cs_vec_norm[2]

p_2_a = dom_1_vec_norm*np.array((0.1/np.sqrt(2),0.05, 0.05))
circle_2_a = np.zeros((3, len(t)))
circle_2_a[0] = p_2_a[0] + r*cos_t*dom_2_vec_norm[0]+r*sin_t*cs_vec_norm[0]
circle_2_a[1] = p_2_a[1] + r*cos_t*dom_2_vec_norm[1]+r*sin_t*cs_vec_norm[1]
circle_2_a[2] = p_2_a[2] + r*cos_t*dom_2_vec_norm[2]+r*sin_t*cs_vec_norm[2]

p_2_b = dom_1_vec_norm*-np.array((0.1/np.sqrt(2),0.05, 0.05))
circle_2_b = np.zeros((3, len(t)))
circle_2_b[0] = p_2_b[0] + r*cos_t*dom_2_vec_norm[0]+r*sin_t*cs_vec_norm[0]
circle_2_b[1] = p_2_b[1] + r*cos_t*dom_2_vec_norm[1]+r*sin_t*cs_vec_norm[1]
circle_2_b[2] = p_2_b[2] + r*cos_t*dom_2_vec_norm[2]+r*sin_t*cs_vec_norm[2]

ax.plot(circle_1_a[0], circle_1_a[1], circle_1_a[2], color = 'b')
ax.plot(circle_1_b[0], circle_1_b[1], circle_1_b[2], color = 'b')
ax.plot(circle_2_a[0], circle_2_a[1], circle_2_a[2], color = 'r')
ax.plot(circle_2_b[0], circle_2_b[1], circle_2_b[2], color = 'r')

ax.set_yticklabels([])
ax.set_xticklabels([])
ax.set_zticklabels([])
ax.grid(False)
plt.axis('off')
ax.view_init(elev=-21, azim=-118)

fig.savefig('/opt/local/l4astro/rbbg94/figures/cs_example.png', overwrite = True)
plt.show()
