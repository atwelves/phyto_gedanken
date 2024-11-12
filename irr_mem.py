# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 11:04:42 2022
Substantial modifications Tue 12 November 2022

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
#Kz = 0.01 # Fast mixing
#Kz = 0.000001 # Slow mixing
# light attenuation coefficient (1/m)
k = 0.04
# mixed layer depth (m)
#D = 20
# time scale to acclimate to irradiance (s)
#gamma = 86400 # Slow acclimation
#gamma = 3600 # Fast acclimation
# maximum irradiance (W/m^2)
I_max = 100
# length of daylight
dlt = 60*60*24
# calculate sine function at dlt
dlt_fn = math.sin(dlt*2*math.pi/86400)
# calculate midpoint of light function
I_avg = I_max*dlt_fn/(dlt_fn-1)
#I_avg = I_max
# number of timesteps
niter = 10080
# number of phytoplankton cells
npart = 10

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

def run_model(D,Kz,gamma):
    # calculate independent trajectory for each phytoplankton cell
    for i in range (0,npart):
        # initialise depth, using uniformly distributed random number
        z[i,0] = random.uniform(0,D)
        # loop over number of timesteps
        for t in range (1,niter):
            # calculate surface irradiance
            I0 = I_avg + (I_max - I_avg)#*math.sin(t*2*math.pi/1440)
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
    return z, I, I_mem 

#for D in (10,20,30,40,50):
count=-1
dimless = numpy.zeros((27))
z_fin   = numpy.zeros((npart,27))
I_fin   = numpy.zeros((npart,27))
m_fin   = numpy.zeros((npart,27))
for D in (20,30,40):
    for Kz_exp in (2,3,4):
        Kz = 1/numpy.power(10,Kz_exp)
        for gamma_exp in (12,24,48):
            count          = count+1
            gamma          = gamma_exp*3600
            #gamma   = 3600*numpy.power(2,gamma_exp)
            dimless[count] = D*D/(Kz*gamma)
            dimless[count] = dimless[count]/(dimless[count]+10)
            results        = run_model(D,Kz,gamma)
            z              = results[0]
            I              = results[1]
            I_mem          = results[2]
            z_fin[:,count] = z[:,-1]
            I_fin[:,count] = I[:,-1]
            m_fin[:,count] = I_mem[:,-1]

print(dimless)

# ----------- #

# Create netCDF files

import netCDF4 as nc
from netCDF4 import Dataset

ncfile = nc.Dataset('lagrange_{}_{}.nc'.format(Kz,gamma),'w', format='NETCDF4')

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




ncfile = nc.Dataset('final_state.nc','w', format='NETCDF4')

# particle axis 
particle_dim = ncfile.createDimension('particle',npart)
# time axis
scale_dim = ncfile.createDimension('scale',27)

ncfile.title='Lagrangian model results'
print(ncfile.title)

particle = ncfile.createVariable('particle',numpy.int32,('particle'))
particle.units = ''
particle.long_name = 'particle identifier'

epsilon = ncfile.createVariable('epsilon',numpy.float64,('scale'))
epsilon.units = ''
epsilon.long_name = 'Ratio of mixing timescale to adaptation timescale'

position = ncfile.createVariable('position',numpy.float64,('particle','scale'))
position.units = 'metres'

irradiance = ncfile.createVariable('irradiance',numpy.float64,('particle','scale'))
irradiance.units = 'W/(m^2)'

irr_memory = ncfile.createVariable('irr_memory',numpy.float64,('particle','scale'))
irr_memory.units = 'W/(m^2)'

# fill data fields

particle[:] = range(npart)
epsilon[:] = dimless[:]

position[:,:] = z_fin
irradiance[:,:] = I_fin
irr_memory[:,:] = m_fin

ncfile.close()

print(dimless)



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


            
        
    
    






