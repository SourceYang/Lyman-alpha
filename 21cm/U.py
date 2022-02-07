#!/usr/bin/env python
# coding: utf-8

# In[4]:


import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
import scipy.special as special
import sys

H_0 = 2.184e-18 # current hubble constant in s^-1
H_r = 3.24076e-18 # 100 km/s/Mpc to 1/s
h_0 = H_0 / H_r
omega_m = 0.315 # current matter abundance
omega_lambda = 0.685 # cosmological constant
G = 6.6743e-8 # gravitational constant in cm^3 g^-1 s^-2
f_He = 0.079 # helium fraction (do we use this number or it changes with redshitf?)
f_H = 0.76 # hydrogen fraction
nu_0 = 2.47e15 * 4 / 3 # center frequency of Lyman-alpha emission, in [Hz]
C = 29979245800.0 # speed of light, in [cm s^-1]
h = 6.6261e-27 # planck constant, in [cm^2 g s^-1]
k = 1.380649e-16 # Boltzmann constant, in [g cm^2 s^-2 K^-1]
omega_b_0 = 0.0224 / h_0 / h_0  # current baryon abundance

z = 8.11 # redshift, could be modified
file_gamma_12 = r'Gamma12aveHII_z008.11_HIIfilter1_RHIImax50_200_300Mpc' # direction of gamma_12
file_Ts = r'Ts_z008.11_L_X3.2e+40_alphaX1.0_TvirminX0.0e+00_zetaIon25.00_Pop2_200_300Mpc' # direction of Ts
file_delta_m = r'updated_smoothed_deltax_z008.11_200_300Mpc' # direction of delta_m


# In[5]:


def load_binary_data(filename, dtype=np.float32):
     """
     We assume that the data was written
     with write_binary_data() (little endian).
     """
     f = open(filename, "rb")
     data = f.read()
     f.close()
     _data = np.frombuffer(data, dtype)
     if sys.byteorder == 'big':
       _data = _data.byteswap()
     return _data


# In[7]:


gamma_12_origin = load_binary_data(file_gamma_12, dtype=np.float32) # photoionization rate
gamma_12 = [x * 1e-12 for x in gamma_12_origin] # convert to cgs unit


# In[ ]:


def sigma_HI(dnu):
    sigma_0 = 6.3e-18 # cm^2
    nu = dnu + nu_0
    epsilon = np.sqrt(nu / nu_0 - 1)
    sigma = sigma_0 * (nu_0 / nu)**4 * np.exp(4 - 4 * np.arctan(epsilon) / epsilon) / (1 - np.exp(-2 * np.pi / epsilon))
    return sigma

def F_nu(dnu, T_bb):
    nu = dnu + nu_0
    f_nu = 2 * h * nu**3 / C**2 / (np.exp(h * nu / k / T_bb) - 1)/ h / nu
    return f_nu

def int_F(T):
    lowerlimit = (2.17872e-11 / h) - nu_0
    upperlimit = (4 * 2.17872e-11 / h) - nu_0
    int_f = integrate.quad(F_nu, lowerlimit, upperlimit, args=T)
    return int_f

def H_F(dnu, T_bb):
    sigma_0 = 6.3e-18 # cm^2
    nu = dnu + nu_0
    epsilon = np.sqrt(nu / nu_0 - 1)
    sigma = sigma_0 * (nu_0 / nu)**4 * np.exp(4 - 4 * np.arctan(epsilon) / epsilon) / (1 - np.exp(-2 * np.pi / epsilon))
    ex = np.float128(h * nu / k / T_bb)
    f_nu = 2 * h * nu**3 / C**2 / (np.exp(ex) - 1) / h / nu
    return sigma * f_nu

def int_Fsigma(T):
    lowerlimit = (2.17872e-11 / h) - nu_0
    upperlimit = (4 * 2.17872e-11 / h) - nu_0
    int_F = integrate.quad(H_F, lowerlimit, upperlimit, args=T)
    return int_F

def sigma_avg(T_bb):
    Fsig = int_Fsigma(T_bb)
    F = int_F(T_bb)
    sigma = Fsig[0] / F[0]
    return sigma


# In[ ]:


T_bb = 50000 # may be modified
sigma_H_avg = sigma_avg(T_bb)


# In[ ]:


omega_b = omega_b_0 * (1 + z)**3 # should be modified, depend on redshift
rho_crit = 3 * H_r * H_r / 8 / np.pi / G
m_H = 1.6735575e-24 # hydrogen atom mass in g
delta_m = load_binary_data(file_delta_m, dtype=np.float32) # matter overdensity


# In[ ]:


def nH(delta_m_i):
    return omega_b * f_H * rho_crit / m_H * (1 + delta_m_i)


# In[ ]:


n_H = []
for i in range(len(delta_m)):
    n_Hi = nH(delta_m[i])
    n_H.append(n_Hi)


# In[ ]:


def U(gamma_12_i, sigma_H_i, n_Hi):
    return gamma_12_i / sigma_H_i / (n_Hi * (1 + f_He) + gamma_12_i / sigma_H_i / C)


# In[ ]:


U_front = []
for i in range(len(gamma_12)):
    U_i = U(gamma_12[i], sigma_H_avg, n_H[i])
    U_front.append(U_i)
    
f = open('test.txt', 'w')
for i in range(len(U_front)):
    f_i = str(U_front[i])
    f.write(f_i+"\n")
f.close()

