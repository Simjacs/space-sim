import numpy as np
import seaborn as sns
from runkut import runkut_step
from utils import calculate_distance

m_earth = 5.972 * 10**24
M_0 = 1.989 * 10**30
astronomical_unit = 1.496 * 10**11
v0 = 3.0 * 10**4
G = 6.67408 * 10**(-11)
step = 1.0
mass_coords = (0.0, 0.0)

# starting conditions, x-y plane
# start at "12 o'clock"
x0 = 0.0
y0 = 1.0 * astronomical_unit
vx0 = v0
vy0 = 0.0

r0 = calculate_distance((x0, y0), mass_coords)
print(r0)
# starting conditions
x_starting_conds = np.array([x0, vx0])
y_starting_conds = np.array([y0, vy0])
print("x0:", x_starting_conds)
print("y0:", y_starting_conds)

# first step
x_change = runkut_step(r0, vx0, M=M_0, h=step)
y_change = runkut_step(r0, vy0, M=M_0, h=step)
x_new_conds = x_starting_conds + x_change
y_new_conds = y_starting_conds + y_change
print("x1:", x_new_conds)
print("y1:", y_new_conds)

# second step
r1 = calculate_distance((x_new_conds[0], y_new_conds[0]), mass_coords)
x_change = runkut_step(r1, x_new_conds[1], M=M_0, h=step)
y_change = runkut_step(r1, y_new_conds[1], M=M_0, h=step)
x_new_conds += x_change
y_new_conds += y_change
print("x2:", x_new_conds)
print("y2:", y_new_conds)

x_conds = np.array([x0, vx0])
y_conds = np.array([y0, vy0])
x_list = []
y_list = []
for i in range(5000):
    print(i)
    r1 = calculate_distance((x_conds[0], y_conds[0]), mass_coords)
    x_change = runkut_step(r1, x_conds[1], M=M_0, h=step)
    y_change = runkut_step(r1, y_conds[1], M=M_0, h=step)
    x_conds += x_change
    y_conds += y_change
    print(f"x{i}:", x_conds)
    print(f"y{i}:", y_conds)
    x_list.append(x_conds[0])
    y_list.append(y_conds[0])

plot = sns.scatterplot(x=x_list, y=y_list)
fig = plot.get_figure()
fig.savefig("fig.png")
