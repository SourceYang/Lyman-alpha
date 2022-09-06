import numpy as np
import sys
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import py21cmfast as p21c
import h5py
import EdgeFinding3dWorkingClass
import intensity_tool
from importlib import reload
import time
import pickle
# reload(EdgeFinding3dWorkingClass2)

parsec=3.085677581e16
H_0=67.74e3/(parsec*10**6) # in units of s^-1, Hubble constants now, 67.74 km/s/mpc
Omega_m=0.3089 # Omega_m = 0.3089+-0.0062
G=6.674e-11  #6.674×10−11 m3*kg−1*s−2 ### 4.30091(25)×10−3 pc*M_solar-1*(km/s)^2
solar_m= 1.98847e30 #(1.98847±0.00007)×10^30 kg
mH_kg = 1.67e-27 # hydrogen mass , kg
rho_c=3*H_0**2/8/np.pi/G # kg/m^3
c_m_s = 3.e8 # light speed in units of m/s
lambda_lya_m = 1216.e-10 # rest frame lyman alpha wavelength in units of m
nu_lya_Hz = c_m_s / lambda_lya_m
H0 = 67.74

def DA(z):
	'''
	unit: Mpc
	'''
	zs = np.linspace(0,z,2000,endpoint=True)
	chi = c_m_s/1000/H0*integrate.simps(1./np.sqrt(1-Omega_m+Omega_m*(1+zs)**3),zs)

	return chi

def H(z):
	'''
	unit : km/s/Mpc
	'''
	H0 = 67.74 # Hubble parameter in km/s/Mpc
	return H0*np.sqrt(1-Omega_m+Omega_m*(1+z)**3)



class load_boxes(object):
	def __init__(self, directory, redshift):
		self.directory = directory
		self.redshift = redshift

		self.z_index = self.box_index_finder(self.directory,self.redshift)
		self.ionized_box = self.load_ionized_box(self.directory,self.z_index)
		self.perturbed_box = self.load_perturbed_box(self.directory,self.z_index)

		self.xH_box = self.ionized_box.xH_box
		self.gamma12_box = self.ionized_box.Gamma12_box
		self.deltam_box = self.perturbed_box.density


	def box_index_finder(self, directory,redshift):
		'''
		find box index with redshift closed to input redshift
		'''
		self.ionized_boxes = list(p21c.cache_tools.query_cache(direc =directory,kind = 'IonizedBox',show=False))
		self.perturbed_field_boxes = list(p21c.cache_tools.query_cache(direc =directory,kind = 'PerturbedField',show=False))
		redshift_arr = np.array([box[1].redshift for box in self.ionized_boxes])

		redshift = 8.0
		z_index = np.abs(redshift_arr-redshift).argmin()

		return z_index

	def load_ionized_box(self, directory, z_index):
		ionized_box = p21c.cache_tools.readbox(direc =directory,fname=self.ionized_boxes[z_index][0])
		return ionized_box

	def load_perturbed_box(self, directory, z_index):
		perturbed_box = p21c.cache_tools.readbox(direc =directory,fname=self.perturbed_field_boxes[z_index][0])
		return perturbed_box



class fronts_finder:
	def __init__(self, xH_box, gamma12_box):
		self.edge = EdgeFinding3dWorkingClass.Edgefinder(xH_box,gamma12_box)


		# self.paris,lines = self.get_pair_line(self.edge)

	def get_pair_line(self, edge):
		pairs = edge.pairs() 
		lines = edge.lines(pairs)
		return pairs,lines

	def get_front(self, edge, pairs, lines):
		fronts = np.array(edge.fronts(pairs,lines)) 
		centers = np.array([front.center for front in fronts])
		uni_center,uni_index = np.unique(centers,axis=0,return_index=True)
		fronts = fronts[uni_index]

		return fronts, uni_center

	def get_gamma(self, edge, fronts):
		gamma12_locs= edge.gamma12loc(fronts)
		gamma12_vals = edge.gamma12val(fronts,gamma12_locs)

		return gamma12_locs, gamma12_vals




class get_intensity(fronts_finder):
	def __init__(self, redshift, T_bb, xH_box, gamma12_box,deltam_box):
		fronts_finder.__init__(self, xH_box, gamma12_box)
		self.redshift = redshift
		self.T_bb = T_bb
		self.xH_box = xH_box
		self.gamma12_box = gamma12_box
		self.deltam_box = deltam_box


		self.dimension = self.get_dimension(self.xH_box)

		self.I_box = np.zeros((self.dimension,self.dimension,self.dimension))
		self.pI_box = np.zeros((self.dimension,self.dimension,self.dimension))

	def get_dimension(self, box):
		dim = box.shape[0]

		return dim

	def get_triangle_geometry(self, front):
		'''
		front is an object from EdgeFinder
		'''
		vertex1 = np.array(front.vertex1)
		vertex2 = np.array(front.vertex2)
		vertex3 = np.array(front.vertex3)
		center = np.array(front.center)
		nhat = front.nhat
		A_Mpc2 = front.area 

		return vertex1, vertex2, vertex3, center, nhat, A_Mpc2

	def get_N_split(self, vertex1, vertex2, vertex3):
		kmax_Mpc = 2
		L_Mpc = np.max([np.sqrt(np.sum((vertex1 - vertex2)**2)), np.sqrt(np.sum((vertex1 - vertex3)**2)), np.sqrt(np.sum((vertex3 - vertex2)**2))])
		N = int(kmax_Mpc * L_Mpc / np.pi) + 1

		return N


	def point_classifier(self, vertex1, vertex2, vertex3):
		N = self.get_N_split(vertex1, vertex2, vertex3)
		points = []
		p_edge = []
		p_corner = []
		p_interior = []
		for i in range(N+1):
		    for j in range(N-i+1):
		        p = vertex1 + i / N * (vertex2-vertex1) + j / N * (vertex3-vertex1)
		        points.append(p)
		        line1_cdt = (i != N and i != 0 and j == 0)
		        line2_cdt = (i == 0 and j != N and j != 0)
		        line3_cdt = (j==(N-i) and i != 0 and i != N)
		        corner1_cdt = (i==0 and j==0)
		        corner2_cdt = (i==0 and j==N)
		        corner3_cdt = (i==N and j==0)
		        if line1_cdt or line2_cdt or line3_cdt:
		            p_edge.append(p)
		        elif corner1_cdt or corner2_cdt or corner3_cdt:
		            p_corner.append(p)
		        else:
		            p_interior.append(p)

		p_edge = np.array(p_edge)
		p_corner = np.array(p_corner)
		p_interior = np.array(p_interior)

		return points, p_edge, p_corner, p_interior


	def get_point_intensity(self, intensity, redshift, A_Mpc2, N):
		V_cell_Mpc3 = 1
		return intensity * (1+redshift) * (c_m_s/1000) / H(redshift) / (nu_lya_Hz / (1+redshift)) * 2 * A_Mpc2 / V_cell_Mpc3 / N**2


	def point_int(self, point):
		x, y, z = point
		if int(x) == self.dimension: x =0
		if int(y) == self.dimension: y =0
		if int(z) == self.dimension: z =0

		return int(x), int(y), int(z)

	def point_vel_nH(self, point, gamma12_val, redshift):
		x, y, z = self.point_int(point)

		delta_m = self.deltam_box[x,y,z]
		front_quantity = intensity_tool.front_velocity(redshift, self.T_bb, gamma12_val, delta_m)

		velocity = front_quantity.U()
		nH = front_quantity.nH()

		return delta_m, velocity, nH

	def point_intensity(self, nhat, redshift, velocity, nH):
		theta = np.arccos(nhat[-1])
		intensity = intensity_tool.I_pI(redshift, theta, velocity, nH)
		I = intensity.I()
		pI = intensity.pI()

		return I, pI

# gamma12_locs, gamma12_vals = intensity_box.get_gamma(intensity_box.edge, fronts)

# intensity_box = get_intensity(redshift, T_bb, xH_box, gamma12_box, deltam_box)

# for i in range(len(fronts)):
# 	gamma12_val = gamma12_vals[i]
# 	vertex1, vertex2, vertex3, center, nhat, A_Mpc2 = intensity_box.get_triangle_geometry(front)
# 	N = intensity_box.get_N_split(vertex1, vertex2, vertex3)
# 	points, p_edge, p_corner, p_interior = intensity_box.point_classifier(vertex1, vertex2, vertex3)

	# for point in points:
	# 	delta_m, velocity, nH = intensity_box.point_vel_nH(point, gamma12_val, redshift)
	# 	I, pI = intensity_box.point_intensity(nhat, redshift, velocity, nH)
	# 	I = intensity_box.get_point_intensity(I, redshift, A_Mpc2, N)
	# 	pI = intensity_box.get_point_intensity(pI, redshift, A_Mpc2, N)
	# 	x, y, z = intensity_box.point_int(point)

        # if np.all(p_corner==point, axis=1).any(): 
        #     intensity_box.I_box[x,y,z] += 1./6 * I
        #     intensity_box.pI_box[x,y,z] += 1./6 * pI


        # elif (p_edge.size!=0 and np.all(p_edge==point, axis=1).any()):
        #     intensity_box.I_box[x,y,z] += 1./2 * I
        #     intensity_box.pI_box[x,y,z] += 1./2 * pI

        # elif (p_interior.size!=0 and np.all(p_interior==point, axis=1).any()):
        #     intensity_box.I_box[x,y,z] += I
        #     intensity_box.pI_box[x,y,z] += pI




















