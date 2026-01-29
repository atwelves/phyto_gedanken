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
dt = 600
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
dlt_fn = math.cos(dlt*2*math.pi/86400)
# calculate midpoint of light function
I_avg = 100#*dlt_fn/(dlt_fn-1)
#I_avg = I_max
# number of timesteps
niter = 1440
# number of phytoplankton cells
npart = 100
# ratio of mixing to adaptation timescale at which
# the irradiance memory is equally influenced by
# instantaneous and depth averaged irradiance
K_scale = 10

# -------------- #

#
#

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
    I_mld = I_avg*(1-numpy.exp(-k*D))/(k*D)
    for i in range (0,npart):
        # initialise depth, using uniformly distributed random number
        z[i,0] = random.uniform(0,D)
        # loop over number of timesteps
        for t in range (1,niter):
            # calculate surface irradiance
            I0 = I_avg + (I_max - I_avg)*math.cos(t*2*math.pi/144)
            if I0 <= 0:
                I0 = 0.00001
            # diffuse to new depth, using normally distributed random number
            delta_z = numpy.random.normal(0,math.sqrt(2*Kz*dt),1)
            delta_z = numpy.nanmax([-D/2,delta_z])
            delta_z = numpy.nanmin([ D/2,delta_z])
            z[i,t] = z[i,t-1] + delta_z
            # reflective boundary condition at sea surface and mixed layer base
            if z[i,t] < 0:
                z[i,t] = - z[i,t]
            elif z[i,t] > D:
                z[i,t] = 2*D - z[i,t] 
            else:
                z[i,t] = z[i,t]
            # calculate instantaneous irradiance
            I[i,t] = I0*numpy.exp(-k*z[i,t])
            #print(I)
            if t == 1:
                # initialise irradiance memory
                I_mem[i,t] = I[i,t]
            else:
                # update irradiance memory using acclimation timescale gamma
                I_mem[i,t] = I_mem[i,t-1] + dt*(I[i,t] - I_mem[i,t-1])/gamma
                # update position memory
                #z_mem[i,t] = - math.log(I_mem[i,t]/I0)/k
    return z, I, I_mem, I_mld 

#for D in (10,20,30,40,50):
count=-1
dimless = numpy.zeros((27))
I_mld_g = numpy.zeros((27))
z_fin   = numpy.zeros((npart,27))
I_fin   = numpy.zeros((npart,27))
m_fin   = numpy.zeros((npart,27))
for D in (10,20,30):
    for Kz_exp in (1,2,3):
        Kz = 1/numpy.power(10,Kz_exp)
        for gamma_hrs in (1,4,12):
            count          = count+1
            gamma          = gamma_hrs*3600
            #gamma   = 3600*numpy.power(2,gamma_exp)
            dimless[count] = D*D/(Kz*gamma)
            dimless_Kz1    = dimless[count]/(dimless[count]+1)
            dimless_Kz10 = dimless[count]/(dimless[count]+10)
            dimless_Kz100 = dimless[count]/(dimless[count]+100)
            results        = run_model(D,Kz,gamma)
            z              = results[0]
            I              = results[1]
            I_mem          = results[2]
            I_mld_g[count] = results[3]
            z_fin[:,count] = z[:,-1]
            I_fin[:,count] = I[:,-1]
            m_fin[:,count] = I_mem[:,-1]
            ### ----------------------------- ###
            plt.figure(figsize=(10,10))
            plt.scatter(I_fin[:,count],m_fin[:,count],s=100,alpha=0.5,color='k')
            #plt.plot(I_fin[:,count] , I_mld_g[count] + dimless_Kz1  *(I_fin[:,count]-I_mld_g[count]),linewidth=2,color=( 27/255,158/255,119/255),label='param, K_zeta=1')
            m, b = numpy.polyfit(I_fin[:,count],m_fin[:,count],1)
            plt.plot(I_fin[:,count], m*I_fin[:,count] + b, label='lagrangian',linewidth=5,color='k')
            plt.plot([0,100] , I_mld_g[count] + dimless_Kz1  *([0,100]-I_mld_g[count]),linestyle='--',linewidth=5,color=( 27/255,158/255,119/255),label='param, K_zeta=1')
            plt.plot([0,100] , I_mld_g[count] + dimless_Kz10 *([0,100]-I_mld_g[count]),linestyle='--',linewidth=5,color=(217/255, 95/255,  2/255),label='param, K_zeta=10')
            plt.plot([0,100] , I_mld_g[count] + dimless_Kz100*([0,100]-I_mld_g[count]),linestyle='--',linewidth=5,color=(117/255,112/255,179/255),label='param, K_zeta=100')
            #plt.title(dimless[count])
            plt.xlim(40,100)
            plt.ylim(40,100)
            plt.xticks([50,75,100],fontsize=40)
            plt.yticks([50,75,100],fontsize=40)
            plt.grid()
            plt.legend(fontsize=30)
            plt.savefig('D{}_Kz-{}_tacc_{}.png'.format(D,Kz_exp,gamma_hrs))


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

irr_mld = ncfile.createVariable('irr_mld',numpy.float64,('scale'))
irr_mld.units = 'W/(m²)'
irr_mld.long_name = 'Average light in mixed layer'

position = ncfile.createVariable('position',numpy.float64,('particle','scale'))
position.units = 'metres'

irradiance = ncfile.createVariable('irradiance',numpy.float64,('particle','scale'))
irradiance.units = 'W/(m^2)'

irr_memory = ncfile.createVariable('irr_memory',numpy.float64,('particle','scale'))
irr_memory.units = 'W/(m^2)'

# fill data fields

particle[:] = range(npart)
epsilon[:] = dimless[:]
irr_mld[:] = I_mld_g[:]
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


            
        
    
    






