import numpy as np

# on earth or not... specify earth radius and rotation
r = 6.371e6
omega = 7.2925e-5

# specify periods (use Inf for a zero frequency, say, Periods = np.Inf)
day = 86400
Period = np.Inf * day
frequency = 2.0 * np.pi / Period

# specify initial k wave number
k_wavenumbers = 4

# time step: make sure results are robust to halving the time step. 1 hour = 3600 s
dt = 3600 

# integration duration (in hours):
integration_time = 24 * 15
Nsteps = round(integration_time*3600/dt)

# Specify initial ray locations
frcy = 45.
frcx = 210.

frcx, frcy = frcx*np.pi/180.0, frcy*np.pi/180.0
frcxm, frcym = frcx*r, r*np.log((1+np.sin(frcy))/np.cos(frcy))
