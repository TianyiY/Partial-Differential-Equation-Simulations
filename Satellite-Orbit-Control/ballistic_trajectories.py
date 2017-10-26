import math
import numpy
import matplotlib.pyplot

h = 0.1 # s
g = 9.81 # m / s2
acceleration = numpy.array([0., -g])
initial_speed = 20. # m / s

def trajectory():
    angles = numpy.linspace(20., 70., 6)

    num_steps = 30
    x = numpy.zeros([num_steps + 1, 2])
    v = numpy.zeros([num_steps + 1, 2])

    for angle in angles:
        # turn degree angle to radius angle
        angle_radius=math.pi/180.*angle
        x[0, 0]=0.     # horizontal position
        x[0, 1]=0.     # vertical position
        v[0, 0]=initial_speed*math.cos(angle_radius)   # horizontal velocity
        v[0, 1]=initial_speed*math.sin(angle_radius)   # vertical velocity

        # Euler forward
        for step in range(num_steps):
            x[step+1]=x[step]+h*v[step]
            v[step+1]=v[step]+h*acceleration  # a is constant
        matplotlib.pyplot.plot(x[:, 0], x[:, 1])

    matplotlib.pyplot.axis('equal')
    axes = matplotlib.pyplot.gca()
    axes.set_xlabel('Horizontal position in m')
    axes.set_ylabel('Vertical position in m')
    matplotlib.pyplot.show()
    return x, v

trajectory()
