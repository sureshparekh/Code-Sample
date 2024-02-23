from vpython import *
#GlowScript 3.1 VPython


G = 6.7e-11
planet = sphere(pos = vector(2e8, 4e8, 0), radius = 5e7, color = color.orange, make_trail = True)
planet.mass = 1e26
planet.p = vector(-1.5e3, 0, 0) * planet.mass
planet.K = mag2(planet.p)/(2*planet.mass)

probe = sphere(pos = vector(-2.0e9, 5.2e8, 0), radius = 1e7, color = color.yellow, make_trail = True)
probe.mass = 1e4
probe.p = vector(5e3, 0,0) * probe.mass
probe.K = mag2(probe.p)/(2*probe.mass)
print("Initial Kinetic Energy of Probe; ", probe.K, ". Initial KE of Planet", planet.K)

dt = 2e2
r = planet.pos - probe.pos
r0 = r

while mag(r) <= mag(r0):
    rate(300)
    r = planet.pos - probe.pos
    F = G * planet.mass * probe.mass * r.hat / mag2(r)
    
    planet.p = planet.p - F*dt
    probe.p = probe.p + F*dt
    probe.K = mag2(probe.p)/(2*probe.mass)
    
    planet.pos = planet.pos + (planet.p / planet.mass) * dt
    probe.pos = probe.pos + (probe.p / probe.mass) * dt
print("Final KE of Probe: ", probe.K, ". Final KE of Planet: ", planet.K)
    
