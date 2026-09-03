import numpy as np
import os
import time as tp

# LICENSE FILE is included in the upper folder
# This script sets up the inputs for the 3D quasi-magneto-static solver in Fortran
# The pre-compiled Fortran executable will be run via this script as well
# SI units unless otherwise stated
# IMPORTANT - field ouputs are cell centered, not a surfaces, the YEE cell positions are solved for and then averaged to cubic centers.

# select the solver and name the inputs txt file to create - both used at bottom of script to submit via command line for convenience
solver='QMS_3D'
inputs_name='inputs.txt'

# generate needed text file. Creates one if it doesn't exist and replaces if it does.
f=open('{}'.format(inputs_name), 'w')

x_size=70
y_size=70
z_size=70

#basic setup info
f.write('{},{},{}\n'.format(x_size,y_size,z_size))     # grid sizes (x,y,z) - in cubes
f.write('0.01,0.01,0.01\n')                            # grid step sizes (x,y,z) - in meters
f.write('0.005\n')                                     #  delta time step - in seconds - no good rule of thumb here for time step relationship to grid step sizes
f.write('0\n')                                         # number of time steps
f.write('1\n')                                         # every this number time steps will be saved for output fields - 1 is default, reduce output fields with higher values
f.write('z,35\n')                                      # 1D slice for output - normal direction and cube height - starts at 1, unlike array below that starts at 0.
f.write('yes\n')                                       # Bulk fields output for entire volume - 'yes' or 'no'
f.write('0,40E-6,0\n')                                 # background static fields in Tesla

#gradient fields - must be traceless and symmetric to be valid inputs (sum diagonal=0 and gxy=gyx for example) - unfortunate restriction
#units of T/m
f.write('0,0.05E-9,0\n')                               # linear gradients for x,y,z directios for Bx
f.write('0.05E-9,0,0\n')                               # linear gradients for x,y,z directios for By
f.write('0,0,0\n')                                     # linear gradients for x,y,z directios for Bz

#Trajectory type 1: v=v_0+0.5*a_0*t with pw on a_0
#This is the motion of the fields - object is relative to that (opposite)
f.write('1\n')                                         # select trajectory x type for built in functions - only '1' supported right now
# if 1 then:
f.write('5,0,0\n')                                     # initial velocity, constant acceleration, number of time steps to accelerate

f.write('1\n')                                         # select trajectory y type for built in functions - only '1' supported right now
# if 1 then:
f.write('0,0,0\n')                                     # initial velocity, constant acceleration, number of time steps to accelerate

f.write('1\n')                                         # select trajectory z type for built in functions - only '1' supported right now
# if 1 then:
f.write('0,0,0\n')                                     # initial velocity, constant acceleration, number of time steps to accelerate

f.write('2\n')                                         # number of materials, don't use material ID zero, it's ignored by the program as a filter
# for each material submit:
f.write('5\n')                                         # material ID number - don't use 0, it defaults to vacuum
f.write('1,1,1\n')                                     # material relative permeability in x,y,z directions
f.write('5E7,5E7,5E7\n')                               # material electrical conductivity (S/m) in x,y,z directions
# for each material submit:
f.write('6\n')                                         # material ID number - don't use 0, it defaults to vacuum
f.write('1,1,1\n')                                     # material relative permeability in x,y,z directions
f.write('0,0,0\n')                                     # material electrical conductivity (S/m) in x,y,z directions

# Build the geometry ################
# geometry file name - will create it below
geom_filename='geom.bin'
# tell the solver the name of the file to look for
f.write('{}\n'.format(geom_filename))                   
# create empty array - zero's will be considered as vacuum and ignored
data=np.zeros((x_size,y_size,z_size), dtype=np.float64)   
# add the user geometry ID number as desired -  recall first cube is at zero in python arrays and last number is unused
#data[12:15,12:15,13]=5
#data[13,13,13]=6
data[20:50,20:50,20:50]=5
data[25:45,25:45,25:45]=6
#data[25:45,10:55,25:45]=6
#data[42:45,25:45,20:50]=5
# save the geometry file in the right format   
data.flatten(order='F').tofile('{}'.format(geom_filename))
#####################################

f.close()

# runs program and clocks the time to run in minutes
start_time=tp.time()
os.system('./{} {}'.format(solver,inputs_name))
end_time=tp.time()
print('Total simulation time is ',(end_time-start_time)/60, ' minutes')