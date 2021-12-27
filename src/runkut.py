def runkut(r0, v0, G, M, h) -> list:
    vk1 = (G * M) / (r0 ** 2)
    rk1 = r0

    vk2 = (G * M) / (r0 + h * (vk1/2))
    rk2 = r0 + h * (rk1/2)

    vk3 = (G * M) / (r0 + h * (vk2/2))
    rk3 = r0 + h * (rk2/2)

    vk4 = (G * M) / (r0 + h * vk3)
    rk4 = r0 + h * rk3

    v1 = (1/6) * h * (vk1 + vk2 + vk3 + vk4)
    r1 = (1/6) * h * (rk1 + rk2 + rk3 + rk4)

    return [r1, v1]


"""
    k1 = f(x, t)
    k2 = f(x + h * (k1/2), t + (h/2))
    k3 = f(x + h * (k2/2), t + (h/2))
    k4 = f(x + h * k3, t + h)
    
    a = GM/r**2
    
    (1/6) * h * (k1 + k2 + k3 + k4)
"""