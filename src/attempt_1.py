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
import copy


def acceleration(r, m2):
    """
    dv/dt depends on position and interacting mass, not time
    :param r: vector between mass being considered and interacting mass
    :param m2: size of interacting mass
    """
    r_hat = r / np.linalg.norm(r)
    dvdt = -(G * m2 / np.linalg.norm(r) ** 2) * r_hat

    return np.nan_to_num(dvdt)


## advance one timestep:
def runkut(rt0, vt0, dt, m):
    """
    implement 4th order runge-kutta for position and velocity
    :param rt0: position vector at start of timestep
    :param vt0: velocity vector at start of timestep
    :param dt: length of timestep
    :return: position and velocity vector at t = t0 + dt
    """
    vk1 = acceleration(rt0, m)
    rk1 = vt0

    vk2 = acceleration(rt0 + (rk1 * (dt / 2)), m)
    rk2 = vt0 + vk1 * (dt / 2)

    vk3 = acceleration(rt0 + (rk2 * (dt / 2)), m)
    rk3 = vt0 + vk2 * (dt / 2)

    vk4 = acceleration(rt0 + (rk3 * dt), m)
    rk4 = vt0 + vk3 * dt

    vdt = vt0 + (1 / 6) * (vk1 + 2 * vk2 + 2 * vk3 + vk4) * dt  # velocity at time t0 + dt
    rdt = rt0 + (1 / 6) * (rk1 + 2 * rk2 + 2 * rk3 + rk4) * dt  # position at time t0 + dt

    return rdt, vdt


# %%
## constants
M0 = 2 * 10 ** 30  ##solar mass
m_e = 6 * 10 ** 24  ## earth mass
R = 1.496 * 10 ** 8  ## astronomical unit
G = 6.67 * 10 ** (-11)
V = 2.978 * 10 ** 4  ## earth avg orbital speed
R_peri = 0.98329 * R
V_peri = 3.029 * 10 ** 4
R_aphe = 1.01671 * R
V_aphe = 2.929 * 10 ** 4

#mass_names = ["sun", "sun2", "earth"]
#masses = [M0, 2*M0, m_e]
mass_names = ["sun", "earth"]
masses = [M0, m_e]
n = len(masses)

## starting condtions:
#r0 = np.array([[-1*R, 0, 0], [1*R, 0, 0], [0, 5*R, 0]])
#v0 = np.array([[0, -10*V, 0], [0, 10*V, 0], [-30*V, 0, 0]])
r0 = np.array([[0, 0, 0], [R, 0, 0]])
v0 = np.array([[0, 0, 0], [0, 32*V, 0]])
r0_peri = np.array([[0, 0, 0], [R_peri, 0, 0]])
v0_peri = np.array([[0, 0, 0], [0, V_peri, 0]])

dt = 0.0001
rt0 = copy.deepcopy(r0_peri)
vt0 = copy.deepcopy(v0_peri)
timesteps = 400

r_data = np.zeros([timesteps + 1, n, 3])
r_data[0] = r0
for i in range(timesteps):  # number of timesteps
    for m_i in range(n):    # n is number of masses
        r_interactions = np.zeros((n, 3))   # 3 indicates dimensionality (3 spatial dims)
        v_interactions = np.zeros((n, 3))
        for m in range(n):
            if m != m_i:
                dist = rt0[m] - rt0[m_i]
                delta_r_i, v_i = runkut(dist, vt0[m], dt, masses[m])
                #print(r_i)   # , v_i)
                r_interactions[m] = rt0[m] + delta_r_i
                v_interactions[m] = v_i

                if i % 100 == 0:
                    print(i, m_i, m)
                    print("r0", rt0[m])
                    print("r:", r_interactions[m])
                    print("v:", v_interactions[m])

        rt0[m] = sum(r_interactions)
        vt0[m] = sum(v_interactions)
    r_data[i + 1] = rt0

m1_x, m1_y = [r_data[i][0, 0] for i in range(timesteps)], [r_data[i][0, 1] for i in range(timesteps)]
m2_x, m2_y = [r_data[i][1, 0] for i in range(timesteps)], [r_data[i][1, 1] for i in range(timesteps)]
#m3_x, m3_y = [r_data[i][2, 0] for i in range(timesteps)], [r_data[i][2, 1] for i in range(timesteps)]

fig = plt.figure()
plt.plot(m1_x, m1_y, markersize=10)
plt.plot(m2_x, m2_y, markersize=2)
fig.savefig(f"orbit_dt_{dt}_steps_{timesteps}.png")
#plt.plot(m3_x, m3_y)



"""
dt = 1
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
    i += 1

plotting_vals = np.array(r_vec)

x = plotting_vals[:, 0]
y = plotting_vals[:, 1]
ax = plt.subplot(1, 1, 1)
ax.plot(x, y)
ax.plot(0, 0, "or")
plt.show()
"""