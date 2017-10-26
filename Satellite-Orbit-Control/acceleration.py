import numpy

earth_mass = 5.97e24 # kg
moon_mass = 7.35e22 # kg
gravitational_constant = 6.67e-11 # N m2 / kg2

# The origin, or (0,0), is at the center of the earth
def acceleration(moon_position, spaceship_position):
    vector_to_moon=moon_position-spaceship_position  # xm
    vector_to_earth=-spaceship_position # xs
    satellite_acceleration=gravitational_constant*(earth_mass*vector_to_earth/numpy.linalg.norm(vector_to_earth)**3+moon_mass*vector_to_moon/numpy.linalg.norm(vector_to_moon)**3)
    return satellite_acceleration