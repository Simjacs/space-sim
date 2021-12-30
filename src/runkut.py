import numpy as np


def grav_acc(r, M=None):
    G = 6.67408 * 10**(-11)
    if M is None:
        M = 1.989 * 10**30
    return (G * M) / (r ** 2)


def runkut_step(r0, v0, M, h) -> np.array:
    vk1 = grav_acc(r0, M)
    rk1 = v0

    vk2 = grav_acc(r0 + h * (vk1/2), M)
    rk2 = v0 + h * (rk1/2)

    vk3 = grav_acc(r0 + h * (vk2/2), M)
    rk3 = v0 + h * (rk2/2)

    vk4 = grav_acc(r0 + h * vk3, M)
    rk4 = v0 + h * rk3

    v1 = (1/6) * h * (vk1 + vk2 + vk3 + vk4)
    r1 = (1/6) * h * (rk1 + rk2 + rk3 + rk4)

    return np.array([r1, v1])


"""
    k1 = f(x, t)
    k2 = f(x + h * (k1/2), t + (h/2))
    k3 = f(x + h * (k2/2), t + (h/2))
    k4 = f(x + h * k3, t + h)
    
    a = GM/r**2
    
    (1/6) * h * (k1 + k2 + k3 + k4)
"""