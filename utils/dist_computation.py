""" A library for distance computations between simple geometric primitives """

import casadi as cs
from tasho.utils import geometry
import numpy as np

def dist_spheres(sphere1, sphere2):

	""" A function to compute the signed-distance between two spheres """

	dist = cs.sumsqr(sphere1['center'] - sphere2['center'])**0.5 - (sphere1['radius'] + sphere2['radius'])

	# dist_vector = dist/cs.norm_1(dist)

	return dist


def dist_sphere_box(sphere, box, vector = False):

	""" A function to compute the distance between a sphere and a box """

	#convert the sphere center to the box coordinates
	sphere_box_coord = geometry.inv_T_matrix(box["tf"])@cs.vertcat(sphere["center"], cs.DM.ones(1))
	diff = sphere_box_coord[0:3]

	# diff = box["tf"][0:3,3] - sphere["center"]

	mins = diff - (box['xyz_len'] + sphere['radius'])
	maxs = -(box['xyz_len'] + sphere['radius']) - diff
	dist_surfaces = cs.vertcat(mins, maxs)
	dist = -100 #a ridiculously small number
	for i in range(6):
		dist = cs.fmax(dist_surfaces[i], dist)

	if not vector:
		return dist
	else:
		for i in range(6):
			if dist_surfaces[i] == dist:
				index = i
				if i == 0:
					dist_vector = cs.vcat([-1, 0, 0])
				elif i == 1:
					dist_vector = cs.vcat([0, -1, 0])
				elif i == 2:
					dist_vector = cs.vcat([0, 0, -1])
				elif i == 3:
					dist_vector = cs.vcat([1, 0, 0])
				elif i == 4:
					dist_vector = cs.vcat([0, 1, 0])
				else:
					dist_vector = cs.vcat([0, 0, 1])

				return dist, dist_vector


def dist_line_segment(line_seg1, line_seg2):

	"""
	A function to compute the distance between a line_segments including sphere
	[Argument]  "line_seg1": body envelope, "line_seg2": others(ex. obstacle, body)
	Component : 'center'(sphere), 'A','B'(line segment), 'type', 'radius'
	"""

	# type assertion
	assert 'type' in line_seg1 and 'type' in line_seg2
	assert line_seg1['type'] == 'line' or line_seg1['type'] == 'sphere' or line_seg1['type'] == 'line_mov' or line_seg1['type'] == 'sphere_mov'
	assert line_seg2['type'] == 'line' or line_seg2['type'] == 'sphere' or line_seg2['type'] == 'line_mov' or line_seg2['type'] == 'sphere_mov'

	# type check -> line or sphere
	if line_seg1['type']=='line' and line_seg2['type']=='line':
		# Between two lines
		pnt_A = line_seg1["A"]
		pnt_B = line_seg1["B"]
		pnt_C = cs.MX(line_seg2["A"])
		pnt_D = cs.MX(line_seg2["B"])
		stp=2
	elif line_seg1['type'] == 'line' and line_seg2['type']=='sphere':
		# line body and sphere obs
		pnt_A = line_seg1["A"]
		pnt_B = line_seg1["B"]
		pnt_C = cs.MX(line_seg2["center"])
		pnt_D = cs.MX(line_seg2["center"])
		anl_u = 0
		stp = 4

	elif line_seg1['type'] == 'line' and line_seg2['type']=='line_mov':
		# line body and moving line obs
		pnt_A = line_seg1["A"]
		pnt_B = line_seg1["B"]
		pnt_C = line_seg2["A"]
		pnt_D = line_seg2["B"]
		stp = 2

	elif line_seg1['type'] == 'line' and line_seg2['type']=='sphere_mov':
		# line body and moving sphere obs
		pnt_A = line_seg1["A"]
		pnt_B = line_seg1["B"]
		pnt_C = line_seg2["center"]
		pnt_D = line_seg2["center"]
		anl_u = 0
		stp = 4

	elif line_seg1['type'] == 'sphere' and line_seg2['type'] == 'line':
		# sphere body and line obs
		pnt_A = cs.MX(line_seg2["A"])
		pnt_B = cs.MX(line_seg2["B"])
		pnt_C = line_seg1["center"]
		pnt_D = line_seg1["center"]
		anl_u = 0
		stp = 4

	elif line_seg1['type'] == 'sphere' and line_seg2['type']=='sphere':
		# sphere body and sphere obs
		pnt_A = line_seg1["center"]
		pnt_B = line_seg1["center"]
		pnt_C = cs.MX(line_seg2["center"])
		pnt_D = cs.MX(line_seg2["center"])
		anl_t = 0
		anl_u = 0
		stp = 5

	elif line_seg1['type'] == 'sphere' and line_seg2['type'] == 'line_mov':
		# sphere body and moving line obs
		pnt_A = line_seg2["A"]
		pnt_B = line_seg2["B"]
		pnt_C = line_seg1["center"]
		pnt_D = line_seg1["center"]
		anl_u = 0
		stp = 4

	elif line_seg1['type'] == 'sphere' and line_seg2['type']=='sphere_mov':
		# sphere body and moving sphere obs
		pnt_A = line_seg1["center"]
		pnt_B = line_seg1["center"]
		pnt_C = line_seg2["center"]
		pnt_D = line_seg2["center"]
		anl_t = 0
		anl_u = 0
		stp = 5

	di1 = pnt_B-pnt_A
	di2 = pnt_D-pnt_C
	di12 = pnt_C-pnt_A
	anl_R = di1.T @ di2
	anl_S1 = di1.T @ di12; anl_S2=di2.T @ di12
	anl_D1 = di1.T @ di1; anl_D2 = di2.T @ di2

	spcs = anl_D1*anl_D2-anl_R**2

	if stp==2:
		# When two lines are parallel, tmp_t is gonna be NaN, because spcs should be zero
		# But, in min-max projection, tmp_t is projected as 0
		tmp_t=(anl_S1*anl_D2-anl_S2*anl_R)/spcs
		anl_t = cs.fmin(1,cs.fmax(0,tmp_t))
		stp+=1
	if stp==3:
		tmp_u=(anl_t*anl_R-anl_S2)/anl_D2
		anl_u = cs.fmin(1,cs.fmax(0,tmp_u))
		stp+=1
	if stp==4:
		tmp_t=(anl_u*anl_R+anl_S1)/anl_D1
		anl_t = cs.fmin(1,cs.fmax(0,tmp_t))
		stp+=1
	if stp==5:
		dist=cs.sqrt(cs.sumsqr(di1*anl_t-di2*anl_u-di12))-(line_seg1['radius']+line_seg2['radius'])

	return dist

if __name__ == '__main__':

	print("no syntax errors")

	#Adding some tests
#	cube = {}
#	cube['tf'] = np.array([[1, 0, 0, 0.5], [0, 1, 0, 0], [0, 0, 1, 0.15], [0, 0, 0, 1]])
#	cube['xyz_len'] = np.array([0.15, 0.15, 0.15])
#	ball = {'center': np.array([0.5, 0.35, 0.15]), 'radius': 0.1}
#	assert cs.fabs(dist_sphere_box(ball, cube) - 0.1) <= 1e-12

#	ball2 = {'center': np.array([0.5, 0.25, 0.15]), 'radius': 0.1}
#	assert cs.fabs(dist_spheres(ball, ball2) + 0.1) <= 1e-12
#	body1={"A":[1,1,2],"B":[2,2,3]}
#	body2={"A":[23,14,12],"B":[20,19,20]}
	#body1={"A":np.array([2,2,1]),"B":np.array([3,3,1])}
	#body2={"A":np.array([1,1,1]),"B":np.array([0,0,1])}
	body1={"A":cs.MX([2,2,1]),"B":cs.MX([3,3,1])}
	body2={"A":[1,1,1],'B':[1,1,10],'obs':True}
	print(dist_line_segment(body1,body2))
