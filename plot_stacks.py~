from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

ang_rad = np.linspace(1,89,89)
no_stacks = 6

for i in range(no_stacks):
	x_corr = np.genfromtxt('/opt/local/l4astro/rbbg94/data/cs_dir_stack_'+str(i)+'.csv', dtype=complex, delimiter = ',')


	cs_lon = 207.8
	cs_lat = -56.3
	cs_vec = hp.pixelfunc.ang2vec(cs_lon, cs_lat, lonlat=True)
	rot_vec = cs_vec/np.linalg.norm(cs_vec)

	vec_lon = 207.8
	vec_lat = cs_lat + 90
	vec = hp.pixelfunc.ang2vec(vec_lon, vec_lat, lonlat=True)

	rot_ang = np.pi/(2*no_stacks)

	dir_1_vec = vec*np.cos(i*rot_ang)+np.cross(rot_vec,vec)*np.sin(i*rot_ang)+rot_vec*np.dot(rot_vec,vec)*(1-np.cos(i*rot_ang))
	dir_1 = hp.pixelfunc.vec2ang(dir_1_vec, lonlat=True)
	dir_2_vec = vec*np.cos(i*rot_ang+np.pi/2)+np.cross(rot_vec,vec)*np.sin(i*rot_ang+np.pi/2)+rot_vec*np.dot(rot_vec,vec)*(1-np.cos(i*rot_ang+np.pi/2))
	dir_2 = hp.pixelfunc.vec2ang(dir_2_vec, lonlat=True)


	fig, ax = plt.subplots(figsize=(18,12))
	ax.axhline(y=0, color = 'black')
	ax.plot(ang_rad, x_corr, color = 'black')
	plt.xticks(np.arange(0, 91, 10))
	plt.xlim(0,90)
	ax.set_xlabel(r'$\alpha/^\circ$')
	ax.set_ylabel('$S$')
	ax.annotate('Directions Used: (lon='+str(round(dir_1[0],1))+'$^\circ$, lat='+str(round(dir_1[1],1))+'$^\circ$), (lon='+str(round(dir_2[0],1))+'$^\circ$, lat='+str(round(dir_2[1],1))+'$^\circ$)', xy = (0.6,0.95), xycoords = 'axes fraction', bbox=dict(facecolor='none', edgecolor='black'))

	fig.savefig('/opt/local/l4astro/rbbg94/figures/corr_0_cs_stacked_'+str(i)+'.png', overwrite = True)

plt.show()
