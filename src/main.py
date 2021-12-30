from runkut import runkut_step
from utils import calculate_distance

m_earth = 5.972 * 10**24
M_0 = 1.989 * 10**30
astronomical_unit = 1.496 * 10**11
v0 = 3 * 10**4
G = 6.67408 * 10**(-11)
step = 1
mass_coords = (0, 0)

# starting conditions, x-y plane
# start at "12 o'clock"
x0 = 0
y0 = 1 * astronomical_unit
vx0 = v0
vy0 = 0

r0 = calculate_distance(mass_coords, (x0, y0))
print(r0)
x1, vx1 = runkut_step(r0, vx0, M=M_0, h=step)
y1, vy1 = runkut_step(r0, vy0, M=M_0, h=step)

print("x:", x1, vx1)
print("y:", y1, vy1)
