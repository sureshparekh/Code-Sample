from vpython import *
#GlowScript 3.1 VPython
G = 6.7e-11

sun = sphere(pos = vector(0, 0, 0), radius = 2e10, color = color.orange, make_trail = True)

sun.mass = 2e30

comet = sphere(pos = vector(-1.5e11, 0, 0), radius = 4e10, color = color.white, make_trail = True)
comet.mass = 1.25e13
comet.p = vector(0, -41.89e3, 0) * comet.mass

sun.p = -comet.p

dt = 1e5
t = 0
while True:
    rate(2000)
    r = comet.pos - sun.pos
    F = G * sun.mass * comet.mass * r.hat / mag2(r)
    
    sun.p = sun.p + F*dt
    comet.p = comet.p - F*dt
    sun.pos = sun.pos + (sun.p/sun.mass) * dt
    comet.pos = comet.pos + (comet.p/comet.mass) * dt
    t = t + dt
