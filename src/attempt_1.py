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
# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def v_derivs(r):
    """
    dv/dt only depends on position, not time
    """
    r_hat = r / np.linalg.norm(r)
    dvdt = -(G * M / np.linalg.norm(r)**2) * r_hat

    return np.nan_to_num(dvdt)

## advance one timestep:
def incr_time(rt0, vt0, dt):
    """
    :param rt0: position vector at start of timestep
    :param vt0: velocity vector at start of timestep
    :param dt: length of timestep
    :return: position and velocity vector at t = t0 + dt
    """
    vk1 = v_derivs(rt0)
    rk1 = vt0
    #print("k1", rk1, vk1)

    vk2 = v_derivs(rt0 + (vk1 * (dt / 2)))
    rk2 = vt0 + vk1 * (dt / 2)
    #print("k2", rk2, vk2)

    vk3 = v_derivs(rt0 + (vk2 * (dt / 2)))
    rk3 = vt0 + vk2 * (dt / 2)
    #print("k3", rk3, vk3)

    vk4 = v_derivs(rt0 + (vk3 * dt))
    rk4 = vt0 + vk3 * dt
    #print("k4", rk4, vk4)

    vdt = vt0 + (1 / 6) * (vk1 + 2 * vk2 + 2 * vk3 + vk4) * dt  # velocity at time t0 + dt
    rdt = rt0 + (1 / 6) * (rk1 + 2 * rk2 + 2 * rk3 + rk4) * dt  # position at time t0 + dt

    #print("dt", rdt, vdt)
    return rdt, vdt

# %%
## constants
M = 2 * 10 ** 30  ##solar mass
R = 1.496 * 10 ** 8  ## astronomical unit
G = 6.67 * 10 ** (-11)
V = 3 * 10 ** 4  ## earth avg orbital speed

## initial conditions
r0 = np.array([R, 0, 0])
v0 = np.array([0, -32*V, 0])

dt = 0.01
rt0 = r0
vt0 = v0
r_vec = []
v_vec = []
print(r0, v0)
for i in range(0, 200000):
    if i % 5000 == 0:
        print(i)
    r_i, v_i = incr_time(rt0, vt0, dt)
    r_vec.append(r_i)
    v_vec.append(v_i)
    #print(i, r_i - rt0, v_i - vt0)

    rt0 = r_i
    vt0 = v_i
    #t += dt



    i += 1

plotting_vals = np.array(r_vec)

x = plotting_vals[:, 0]
y = plotting_vals[:, 1]
ax = plt.subplot(1, 1, 1)
ax.plot(x, y)
ax.plot(0,0, "or")
plt.show()
