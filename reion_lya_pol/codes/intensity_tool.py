"""
Adapted from Yuanyuan's interpolation codes
"""


import math
import numpy as np
import scipy

import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.special as special
import sys
import h5py

H_0 = 2.184e-18 # current hubble constant in s^-1
H_r = 3.24076e-18 # 100 km/s/Mpc to 1/s
h_0 = H_0 / H_r
omega_m = 0.315 # current matter abundance
omega_lambda = 0.685 # cosmological constant
G = 6.6743e-8 # gravitational constant in cm^3 g^-1 s^-2
f_He = 0.079 # helium fraction
f_H = 0.76 # hydrogen fraction
nu_0 = 2.47e15 * 4 / 3 # center frequency of Lyman-alpha emission, in [Hz]
C = 29979245800.0 # speed of light, in [cm s^-1]
h = 6.6261e-27 # planck constant, in [cm^2 g s^-1]
k = 1.380649e-16 # Boltzmann constant, in [g cm^2 s^-2 K^-1]
omega_b_0 = 0.0224 / h_0 / h_0  # current baryon abundance
m_H = 1.6735575e-24 # hydrogen atom mass in g

T_bb = 47834.5 # 21cmFAST setting blackbody temperature in K

vlya = 2.47e15 # center frequency of Lyman alpha photon, in [Hz]


"""
To use this class, call I_pI(z, theta, Ui, nHi)

z is redshift, theta is the direction of view from 0 to 2pi, Ui is front velocity in cm/s,
nHi is the number density of neutral hydrogen atom in cm^-3

all variables should be individual values rather than array or list

.I() will return to the intensity
.pI() will return to the polarized intensity
"""

class I_pI:
    def __init__(self, z, theta, Ui, nHi):
        self.z = z
        self.theta = theta
        self.mu = np.cos(theta)
        self.Ui = Ui
        self.nHi = nHi
        self.U = []
        self.nH = []
        self.read_data()

    def read_data(self):
        with open('./pmu_m.txt') as f:
            f = [x.strip() for x in f if x.strip()]
            data = [tuple(map(float,x.split())) for x in f[0:]]
            self.index_U = [x[0] for x in data]
            self.index_nH = [x[1] for x in data]
            self.mean_p = [x[2] for x in data]
            self.std_p = [x[3] for x in data]
            self.mean_n = [x[4] for x in data]
            self.std_n = [x[5] for x in data]

        with open(r'./b.txt') as f:
            f = [x.strip() for x in f if x.strip()]
            data = [tuple(map(float,x.split())) for x in f[0:]]
            index_U = [x[0] for x in data]
            index_nH = [x[1] for x in data]
            self.b0 = [x[2] for x in data]
            self.b1 = [x[3] for x in data]
            self.b2 = [x[4] for x in data]
            self.b3 = [x[5] for x in data]
            self.b4 = [x[6] for x in data]
            self.b5 = [x[7] for x in data]
    
        with open(r'./c.txt') as f:
            f = [x.strip() for x in f if x.strip()]
            data = [tuple(map(float,x.split())) for x in f[0:]]
            index_U = [x[0] for x in data]
            index_nH = [x[1] for x in data]
            self.c2 = [x[2] for x in data]
            self.c3 = [x[3] for x in data]
            self.c4 = [x[4] for x in data]
            self.c5 = [x[5] for x in data]
            self.c6 = [x[6] for x in data]
            self.c7 = [x[7] for x in data]

        for i in range(len(self.index_U)):
            if self.index_U[i] == 11:
                u_i = 7e6
            if self.index_U[i] == 12:
                u_i = 1e7
            if self.index_U[i] == 13:
                u_i = 3e7
            if self.index_U[i] == 14:
                u_i = 7e7
            if self.index_U[i] == 15:
                u_i = 1e8
            if self.index_U[i] == 16:
                u_i = 3e8
            if self.index_U[i] == 17:
                u_i = 7e8
            if self.index_U[i] == 18:
                u_i = 1e9
            if self.index_U[i] == 19:
                u_i = 3e9
            if self.index_U[i] == 20:
                u_i = 7e9
            if self.index_U[i] == 21:
                u_i = 1e10
            if self.index_U[i] == 22:
                u_i = 2.7e10
            self.U.append(float(u_i))
    
            if self.index_nH[i] == 0:
                nH_i = 1e-9
            if self.index_nH[i] == 1:
                nH_i = 5e-9
            if self.index_nH[i] == 2:
                nH_i = 1e-8
            if self.index_nH[i] == 3:
                nH_i = 5e-8
            if self.index_nH[i] == 4:
                nH_i = 1e-7
            if self.index_nH[i] == 5:
                nH_i = 5e-7
            if self.index_nH[i] == 6:
                nH_i = 1e-6
            if self.index_nH[i] == 7:
                nH_i = 5e-6
            if self.index_nH[i] == 8:
                nH_i = 1e-5
            if self.index_nH[i] == 9:
                nH_i = 5e-5
            if self.index_nH[i] == 10:
                nH_i = 1e-4
            if self.index_nH[i] == 11:
                nH_i = 5e-4
            if self.index_nH[i] == 12:
                nH_i = 1e-3
            if self.index_nH[i] == 13:
                nH_i = 5e-3
            if self.index_nH[i] == 14:
                nH_i = 1e-2
            if self.index_nH[i] == 15:
                nH_i = 5e-2
            if self.index_nH[i] == 16:
                nH_i = 1e-1
            if self.index_nH[i] == 17:
                nH_i = 5e-1
            if self.index_nH[i] == 18:
                nH_i = 1
            self.nH.append(float(nH_i))
    
    def fit_2d(self, x, y, z):
        if x < 1e7:
            U_index = 11
        if 1e7 <= x < 3e7:
            U_index = 12
        if 3e7 <= x < 7e7:
            U_index = 13
        if 7e7 <= x < 1e8:
            U_index = 14
        if 1e8 <= x < 3e8:
            U_index = 15
        if 3e8 <= x < 7e8:
            U_index = 16
        if 7e8 <= x < 1e9:
            U_index = 17
        if 1e9 <= x < 3e9:
            U_index = 18
        if 3e9 <= x < 7e9:
            U_index = 19
        if 7e9 <= x <1e10:
            U_index = 20
        if x >= 1e10:
            U_index = 21
        
        if y < 5e-9:
            nH_index = 0
        if 5e-9 <= y < 1e-8:
            nH_index = 1
        if 1e-8 <= y < 5e-8:
            nH_index = 2
        if 5e-8 <= y < 1e-7:
            nH_index = 3
        if 1e-7 <= y < 5e-7:
            nH_index = 4
        if 5e-7 <= y < 1e-6:
            nH_index = 5
        if 1e-6 <= y < 5e-6:
            nH_index = 6
        if 5e-6 <= y < 1e-5:
            nH_index = 7
        if 1e-5 <= y < 5e-5:
            nH_index = 8
        if 5e-5 <= y < 1e-4:
            nH_index = 9
        if 1e-4 <= y < 5e-4:
            nH_index = 10
        if 5e-4 <= y < 1e-3:
            nH_index = 11
        if 1e-3 <= y < 5e-3:
            nH_index = 12
        if 5e-3 <= y < 1e-2:
            nH_index = 13
        if 1e-2 <= y < 5e-2:
            nH_index = 14
        if 5e-2 <= y < 1e-1:
            nH_index = 15
        if 1e-1 <= y < 5e-1:
            nH_index = 16
        if y >= 5e-1:
            nH_index = 17
            
        j = self.index_U.index(U_index)
        j1 = self.index_U.index(U_index+1)
        
        jk = self.index_nH.index(nH_index, j, j+19)
        j1k = self.index_nH.index(nH_index, j1, j1+19)
        jk1 = self.index_nH.index(nH_index+1, j, j+19)
        j1k1 = self.index_nH.index(nH_index+1, j1, j1+19)
    
        f_A = z[jk] + (x - self.U[jk]) / (self.U[j1k] - self.U[jk]) * (z[j1k] - z[jk])
        f_B = z[jk1] + (x - self.U[jk]) / (self.U[j1k] - self.U[jk]) * (z[j1k1] - z[jk1])
        f_xy = f_A + (y - self.nH[jk]) / (self.nH[jk1] - self.nH[jk]) * (f_B - f_A)
    
        return f_xy


    def P0(self, x):
        return 1
    def P1(self, x):
        return x
    def P2(self, x):
        return (3*x**2)/2. - 1/2.
    def P3(self, x):
        return (5*x**3)/2. - 3*x/2.
    def P4(self, x):
        return (35*x**4)/8. - (15*x**2)/4. + 3/8.
    def P5(self, x):
        return (63*x**5)/8. - (35*x**3)/4. + 15*x/8.

    def P22(self, x):
        return 3*(1-x**2)
    def P23(self, x):
        return 15*x*(1-x**2)
    def P24(self, x):
        return 15/2 * (7*x**2-1)*(1-x**2)
    def P25(self, x):
        return 1/4*(11*x*self.P24(x) - 7*self.P23(x))
    def P26(self, x):
        return 1/5*(13*x*self.P25(x) - 8*self.P24(x))
    def P27(self, x):
        return 1/6*(15*x*self.P26(x) - 9*self.P25(x))
    
    def b0_func(self):
        return self.fit_2d(self.Ui, self.nHi, self.b0)

    def b1_func(self):
        return self.fit_2d(self.Ui, self.nHi, self.b1)
    def b2_func(self):
        return self.fit_2d(self.Ui, self.nHi, self.b2)
    def b3_func(self):
        return self.fit_2d(self.Ui, self.nHi, self.b3)
    def b4_func(self):
        return self.fit_2d(self.Ui, self.nHi, self.b4)
    def b5_func(self):
        return self.fit_2d(self.Ui, self.nHi, self.b5)

    def pmu(self):
        b0_fit = self.fit_2d(self.Ui, self.nHi,self.b0)
        b1_fit = self.fit_2d(self.Ui, self.nHi, self.b1)
        b2_fit = self.fit_2d(self.Ui, self.nHi, self.b2)
        b3_fit = self.fit_2d(self.Ui, self.nHi, self.b3)
        b4_fit = self.fit_2d(self.Ui, self.nHi, self.b4)
        b5_fit = self.fit_2d(self.Ui, self.nHi, self.b5)
        y = self.mu
        pmu_fit = b0_fit + b1_fit * self.P1(y) + b2_fit * self.P2(y) + b3_fit * self.P3(y) + b4_fit * self.P4(y) + b5_fit * self.P5(y)
        
        return pmu_fit

    def n_rate(self):
        n_fit = self.fit_2d(self.Ui, self.nHi, self.mean_n)
        return n_fit

    def I(self):
        b0_fit = self.fit_2d(self.Ui, self.nHi, self.b0)
        b1_fit = self.fit_2d(self.Ui, self.nHi, self.b1)
        b2_fit = self.fit_2d(self.Ui, self.nHi, self.b2)
        b3_fit = self.fit_2d(self.Ui, self.nHi, self.b3)
        b4_fit = self.fit_2d(self.Ui, self.nHi, self.b4)
        b5_fit = self.fit_2d(self.Ui, self.nHi, self.b5)
        y = self.mu
        pmu_fit = b0_fit + b1_fit * self.P1(y) + b2_fit * self.P2(y) + b3_fit * self.P3(y) + b4_fit * self.P4(y) + b5_fit * self.P5(y)
        
        n_fit = self.fit_2d(self.Ui, self.nHi, self.mean_n)
        
        vlya_z = vlya / (1 + self.z)
        
        return h * vlya_z * n_fit * pmu_fit / (2 * np.pi  * (1 + self.z)**3)
    
    def pI(self):
        c2_fit = self.fit_2d(self.Ui, self.nHi, self.c2)
        c3_fit = self.fit_2d(self.Ui, self.nHi, self.c3)
        c4_fit = self.fit_2d(self.Ui, self.nHi, self.c4)
        c5_fit = self.fit_2d(self.Ui, self.nHi, self.c5)
        c6_fit = self.fit_2d(self.Ui, self.nHi, self.c6)
        c7_fit = self.fit_2d(self.Ui, self.nHi, self.c7)
        y = self.mu
        pImu_fit = c2_fit * self.P22(y) + c3_fit * self.P23(y) + c4_fit * self.P24(y) + c5_fit * self.P25(y) + c6_fit * self.P26(y) + c7_fit * self.P27(y)
        
        n_fit = self.fit_2d(self.Ui, self.nHi, self.mean_n)
        
        vlya_z = vlya / (1 + self.z)
        
        return h * vlya_z * n_fit * pImu_fit / (2 * np.pi * (1 + self.z)**3)


"""
To use this class, call front_velocity(z, T_bb, gamma_12, delta_m)

z is redshift, T_bb is blackbody temperature in K, gamma_12 is the photoionization rate in unit cgs*1e12,
delta_m is the matter overdensity in cgs unit

all variables should be individual values rather than array or list

.nH() will return to the number density of neutral hydrogen atom in cm^-3
.U() will return to the front velocity in cm/s
"""

class front_velocity:
    def __init__(self, z, T_bb, gamma_12, delta_m):
        self.gamma_12 = gamma_12 * 1e-12 # convert to cgs unit
        self.z = z
        self.T_bb = T_bb
        self.delta_m = delta_m
        
    def sigma_HI(self, dnu):
        sigma_0 = 6.3e-18 # cm^2
        nu = dnu + nu_0
        epsilon = np.sqrt(nu / nu_0 - 1)
        sigma = sigma_0 * (nu_0 / nu)**4 * np.exp(4 - 4 * np.arctan(epsilon) / epsilon) / (1 - np.exp(-2 * np.pi / epsilon))
        return sigma
    
    def F_nu(self, dnu, T_bb):
        nu = dnu + nu_0
        f_nu = 2 * h * nu**3 / C**2 / (np.exp(h * nu / k / T_bb) - 1)/ h / nu
        return f_nu

    def int_F(self, T):
        lowerlimit = (2.17872e-11 / h) - nu_0
        upperlimit = (4 * 2.17872e-11 / h) - nu_0
        int_f = integrate.quad(self.F_nu, lowerlimit, upperlimit, args=T)
        return int_f

    def H_F(self, dnu, T_bb):
        sigma_0 = 6.3e-18 # cm^2
        nu = dnu + nu_0
        epsilon = np.sqrt(nu / nu_0 - 1)
        sigma = sigma_0 * (nu_0 / nu)**4 * np.exp(4 - 4 * np.arctan(epsilon) / epsilon) / (1 - np.exp(-2 * np.pi / epsilon))
        ex = np.float64(h * nu / k / T_bb) # if ex larger than 700, should use float128
        f_nu = 2 * h * nu**3 / C**2 / (np.exp(ex) - 1) / h / nu
        return sigma * f_nu

    def int_Fsigma(self, T):
        lowerlimit = (2.17872e-11 / h) - nu_0 # 2.17872e-11 is 13.6eV in cgs unit
        upperlimit = (4 * 2.17872e-11 / h) - nu_0
        int_F = integrate.quad(self.H_F, lowerlimit, upperlimit, args=T)
        return int_F

    def sigma_avg(self, T):
        Fsig = self.int_Fsigma(T)
        F = self.int_F(T)
        sigma = Fsig[0] / F[0]
        return sigma
        
    def nH(self):
        omega_b = omega_b_0 * (1 + self.z)**3
        rho_crit = 3 * H_0 * H_0 / 8 / np.pi / G
        return omega_b * f_H * rho_crit / m_H * (1 + self.delta_m)
        
    def U(self):
        sigma_H_avg = self.sigma_avg(self.T_bb)
        nH = self.nH()
        return self.gamma_12 / sigma_H_avg / (nH * (1 + f_He) + self.gamma_12 / sigma_H_avg / C)
