from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

CMB_DIST = 14000 #Mpc
CELL_SIZE = 320 #Mpc

def circle_finder(cell_size, x, y):
	d = cell_size*(np.sqrt(x**2+y**2))
	ang_rad = np.arccos(d/(2*CMB_DIST))
	return ang_rad

def draw_circle(x, y, r, coords, angles):
	for i in range(len(angles)):
		coords[0] = x + r*np.cos(angles)
		coords[1] = y + r*np.sin(angles)
	return coords 

if __name__ == '__main__':

	coords_0 = np.zeros((2, 361))
	coords_x1 = np.zeros((2, 361))
	coords_x40 = np.zeros((2, 361))
	coords_y40 = np.zeros((2, 361))
	angles = np.radians(np.arange(0,361, 1))

	circle_0 = draw_circle(0, 0, CMB_DIST, coords_0, angles)
	circle_x1 = draw_circle(1*CELL_SIZE, 0, CMB_DIST, coords_x1, angles)
	circle_x40 = draw_circle(40*CELL_SIZE, 0, CMB_DIST, coords_x40, angles)
	circle_y40 = draw_circle(0, 40*CELL_SIZE, CMB_DIST, coords_y40, angles)

	no_of_circles = int(CMB_DIST/CELL_SIZE)*2
	x=0
	y=0

	ang_rad = np.zeros(no_of_circles)

	for i in range(no_of_circles):
		ang_rad[i] = circle_finder(CELL_SIZE, x, y)
		x+=1


	fig, ax = plt.subplots(figsize = (12,8))
	"""ax.plot(circle_0[0], circle_0[1])
	ax.plot(circle_x1[0], circle_x1[1])
	ax.plot(CMB_DIST*np.cos(ang_rad[1]), CMB_DIST*np.sin(ang_rad[1]), 'x', color = 'r')
	ax.plot(CMB_DIST*np.cos(ang_rad[1]), -CMB_DIST*np.sin(ang_rad[1]), 'x', color = 'r')
	ax.set_title('Angular radius = ' + str(np.degrees(ang_rad[1])) + "$^\circ$")"""
	ax.plot(circle_0[0], circle_0[1])
	ax.plot(circle_x40[0], circle_x40[1])
	ax.plot(CMB_DIST*np.cos(ang_rad[40]), CMB_DIST*np.sin(ang_rad[40]), 'x', color = 'r')
	ax.plot(CMB_DIST*np.cos(ang_rad[40]), -CMB_DIST*np.sin(ang_rad[40]), 'x', color = 'r')
	ax.set_title('Angular radius = ' + str(np.degrees(ang_rad[40])) + "$^\circ$")
	"""ax.plot(circle_0[0], circle_0[1])
	ax.plot(circle_y40[0], circle_y40[1])
	ax.plot(CMB_DIST*np.cos(ang_rad[40]), CMB_DIST*np.sin(ang_rad[40]), 'x', color = 'r')
	ax.plot(CMB_DIST*np.cos(ang_rad[40]), -CMB_DIST*np.sin(ang_rad[40]), 'x', color = 'r')
	ax.set_title('Angular radius = ' + str(np.degrees(ang_rad[40])) + "$^\circ$")"""
	plt.tight_layout()
	plt.show()
