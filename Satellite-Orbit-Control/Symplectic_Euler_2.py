import math
import numpy
import matplotlib.pyplot

h = 0.05  # s
g = 9.81  # m / s2
length = 1.  # m

axes_x = matplotlib.pyplot.subplot(311)
axes_x.set_ylabel('x in m')
axes_v = matplotlib.pyplot.subplot(312)
axes_v.set_ylabel('v in m/s')
axes_v.set_xlabel('t in s')
axes_phase_space = matplotlib.pyplot.subplot(313)
axes_phase_space.set_xlabel('x in m')
axes_phase_space.set_ylabel('v in m/s')

def acceleration(position):
    return -g*math.sin(position/length)

def symplectic_euler():

    num_steps = 80
    x = numpy.zeros(num_steps + 1)  # m around circumference
    v = numpy.zeros(num_steps + 1)  # m / s
    colors = [(0, 'g'), (3, 'c'), (10, 'b'), (30, 'm'), (79, 'r')]
    times = h * numpy.arange(num_steps + 1)

    num_initial_conditions = 50

    for i in range(num_initial_conditions):
        phi = 2.*math.pi*i/num_initial_conditions
        x[0] = 2.0+0.25* math.cos(phi)
        v[0] = 2.0* math.sin(phi)

        for step in range(num_steps):
            x[step+1]=x[step]+h*v[step]
            v[step+1]=v[step]+h*acceleration(x[step+1])

        axes_x.plot(times, x, c='k', alpha=0.1)
        axes_v.plot(times, v, c='k', alpha=0.1)
        for step, color in colors:
            #matplotlib.pyplot.hold(True)
            axes_x.scatter(times[step], x[step], c=color, edgecolors='none')
            axes_v.scatter(times[step], v[step], c=color, edgecolors='none')
            axes_phase_space.scatter(x[step], v[step], c=color, edgecolors='none', s=4)
            matplotlib.pyplot.show()
    return x, v


symplectic_euler()
