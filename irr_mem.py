# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 11:04:42 2022

@author: agt_9
"""

# script to illustrate interplay between photo-adaptation and mixing timescales
# based on model of Kida & Ito (2017) and first published using MATLAB

# started 16/08/2022

import random
import numpy
import math
import matplotlib.pyplot as plt

# -------------- #

# Declare constants

# time step (s)
dt = 60
# diffusivity parameter (m^2/s)
Kz = 0.01 # Fast mixing
#Kz = 0.001 # Slow mixing
# light attenuation coefficient (1/m)
k = 0.04
# mixed layer depth (m)
D = 20
# time scale to acclimate to irradiance (s)
#gamma = 86400 # Slow acclimation
gamma = 3600 # Fast acclimation
# maximum irradiance (W/m^2)
I_max = 100
# length of daylight
dlt = 60*60*0
# calculate sine function at dlt
dlt_fn = math.sin(dlt*2*math.pi/86400)
# calculate midpoint of light function
I_avg = I_max*dlt_fn/(dlt_fn-1)
#I_avg = I_max
# number of timesteps
niter = 10080
# number of phytoplankton cells
npart = 100

# -------------- #

# Declare variables

# phytoplankton cell position
z = numpy.array([[0.0]*niter]*npart)
# position memory
z_mem = numpy.array([[0.0]*niter]*npart)
# instantaneous irradiance
I = numpy.array([[0.0]*niter]*npart)
# irradiance memory
I_mem = numpy.array([[0.0]*niter]*npart)

# -------------- #

# Run model

# calculate independent trajectory for each phytoplankton cell
for i in range (0,npart):
    # initialise depth, using uniformly distributed random number
    z[i,0] = random.uniform(0,D)
    # loop over number of timesteps
    for t in range (1,niter):
        # calculate surface irradiance
        I0 = I_avg + (I_max - I_avg)*math.sin(t*2*math.pi/1440)
        if I0 <= 0:
            I0 = 0.00001
        # diffuse to new depth, using normally distributed random number
        z[i,t] = z[i,t-1] + numpy.random.normal(0,math.sqrt(2*Kz*dt),1)
        # reflective boundary condition at sea surface and mixed layer base
        if z[i,t] < 0:
            z[i,t] = - z[i,t]
        elif z[i,t] > D:
            z[i,t] = 2*D - z[i,t] 
        else:
            z[i,t] = z[i,t]
        # calculate instantaneous irradiance
        I[i,t] = I0*numpy.exp(-k*z[i,t])
        if t == 1:
            # initialise irradiance memory
            I_mem[i,t] = I[i,t]
        else:
            # update irradiance memory using acclimation timescale gamma
            I_mem[i,t] = I_mem[i,t-1] + (I[i,t] - I_mem[i,t-1])/gamma
            # update position memory
            z_mem[i,t] = - math.log(I_mem[i,t]/I0)/k
            
# ----------- #

# Create netCDF files

import netCDF4 as nc
from netCDF4 import Dataset

ncfile = nc.Dataset('lagrange.nc','w', format='NETCDF4')

# particle axis 
particle_dim = ncfile.createDimension('particle',npart)
# time axis
time_dim = ncfile.createDimension('time',niter)

ncfile.title='Lagrangian model results'
print(ncfile.title)

particle = ncfile.createVariable('particle',numpy.int32,('particle'))
particle.units = ''
particle.long_name = 'particle identifier'

time = ncfile.createVariable('time',numpy.int32,('time'))
time.units = 'seconds'
time.long_name = 'time elapsed'

position = ncfile.createVariable('position',numpy.float64,('particle','time'))
position.units = 'metres'

irradiance = ncfile.createVariable('irradiance',numpy.float64,('particle','time'))
irradiance.units = 'W/(m^2)'

irr_memory = ncfile.createVariable('irr_memory',numpy.float64,('particle','time'))
irr_memory.units = 'W/(m^2)'

# fill data fields

particle[:] = range(npart)
time[:] = numpy.linspace(0,dt*niter,niter)

position[:,:] = z
irradiance[:,:] = I
irr_memory[:,:] = I_mem

ncfile.close()

# ----------- #

# Plot graphs

# declare time axis

time_mins = numpy.linspace(0, 1, 10080)
            
fig, ax = plt.subplots()

ax.plot(time_mins, z[0], linewidth=2.0)

plt.show()
            
fig, ax = plt.subplots()

plt.scatter(time_mins, I_mem[0], linewidth=2.0)

plt.show()


            
        
    
    






