## test scenario: simulation for a pendulum under earth gravity
import numpy as np
import matplotlib.pyplot as plt

"""
d2theta/dt2 = -g/L sin(theta)
reducing this second order eqtn to a system of two first order eqtns gives:
d(omega)/dt = mu * omega - g/l sin(theta)
d(theta)/dt = omega
where g = 9.8 ms**-2
theta is angle from the vertical in radians
and L is the length of the pendulum in metres 
"""


def dwdt(theta, L, g):
    return - (g / L) * np.sin(theta)


def runkut(theta_0, omega_0, dt, g, L):
    omega_k1 = dwdt(theta_0, g, L)
    theta_k1 = omega_0

    omega_k2 = dwdt(theta_0 + (theta_k1 * (dt / 2)), g, L)
    theta_k2 = omega_0 + omega_k1 * (dt / 2)

    omega_k3 = dwdt(theta_0 + (theta_k2 * (dt / 2)), g, L)
    theta_k3 = omega_0 + omega_k2 * (dt / 2)

    omega_k4 = dwdt(theta_0 + (theta_k3 * dt), g, L)
    theta_k4 = omega_0 + omega_k3 * dt

    omega_dt = omega_0 + (1 / 6) * (omega_k1 + 2 * omega_k2 + 2 * omega_k3 + omega_k4) * dt  # omega at time t0 + dt
    theta_dt = theta_0 + (1 / 6) * (theta_k1 + 2 * theta_k2 + 2 * theta_k3 + theta_k4) * dt  # theta at time t0 + dt

    return theta_dt, omega_dt


# constants
g = 9.81
L = 0.5
dt = 0.01
timesteps = 2000

# initial conditions:
theta_t0 = np.pi / 6  # 30 degrees
omega_t0 = 0

theta_values = [theta_t0]
omega_values = [omega_t0]

for i in range(timesteps):
    theta0 = theta_values[-1]
    omega0 = omega_values[-1]
    theta_dt, omega_dt = runkut(theta0, omega0, dt, g, L)

    theta_values.append(theta_dt)
    omega_values.append(omega_dt)

fig = plt.figure()
plt.plot(theta_values, omega_values, "b.")
fig.savefig(f"pendulum_timesteps{timesteps}_dt_{dt}.png")
