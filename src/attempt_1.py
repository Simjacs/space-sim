#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 13:33:12 2020

@author: simran
"""

"""
runge-kutta algorithm:
    y(t + dt) = y(t) + (1/6)*(k1 + 2*k2 + 2*k3 + k4)*dt
    
where:
    k1 = f(y(t), t) ##value at start
    k2 = f(y(t) + k1 * (dt/2), t+(dt/2)) ##estimate 1 for value at midpoint
    k3 = f(y(t) + k2 * (dt/2), t+(dt/2)) ##estimate 2 for value at midpoint
    k4 = f(y(t) + k3 * dt, t+dt) ##estimate of value at end of timestep
    
    f(y(t), t) = y'(t)
"""


"""
euler method:
    y(t + dt) = y(t) + f(y, t)*dt
"""

## need to def ode (eg dydx = x + y) and pass to runkut algorithm
## dxdt = v
## dvdt = (-GM/r**2) * r_hat

"""
Orbit problem

2nd order ode: d2r/dt2 = F/m

radial force depends only on position
F = (-GMm/r**2) * r_hat
r_hat = r/|r|
|r| = np.linalg.norm(r)

"""

"""
3d problem:
    r = [x,y,z]
    v = [vx, vy, vz]
    
    set M at origin
"""
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%%
## constants
M = 2 * 10 ** 30 ##solar mass
m = 6 * 10 ** 24 ## earth mass
R = 1.496 * 10 ** 8 ## astronomical unit
G = 6.67 * 10 ** (-11)
V = 3 * 10**4 ## earth avg orbital speed

## initial conditions
r0 = np.array([R, 0, 0])
v0 = np.array([0, V, 0])


def runkut(r0, v0, t0, dt, v_derivs):
    """
    y0 is value of y at start of timestep
    derivs must be funct of y,t 
    """
    ## runge kutta constants:
    vk1 = derivs(v0, t0)
    rk1 = v0
    
    vk2 = derivs(v0 + (vk1 * (dt/2)), t0 + (dt/2))
    rk2 = v0 + (rk1*(dt/2))
    #rk2 = drdt(r0 + (rk1 * (dt/2)), t0 + (dt/2))
    
    vk3 = derivs(v0 + (vk2 * (dt/2), t0 + (dt/2)))
    ## rk2 = drdt(r0 + (rk2 * (dt/2)), t0 + (dt/2))
    rk2 = v0 + 
    
    vk4 = derivs(v0 + (vk3 * dt), t0 + dt)
    
    ## y(t+dt)
    return(r0 + (1/6)*(vk1 + 2*vk2 + 2*vk3 + vk4) * dt)
    
def dvdt(r, t):
    """
    dv/dt only depends on position, not time
    """    
    r_hat = r/np.linalg.norm(r)
    
    dvdt_vec = -(G*M/r**2) * r_hat
    
    return(np.nan_to_num(dvdt_vec))             
    
 
#%% 
    
def drdt(r, t):
    """
    dr/dt is just v
    """

#%% 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

