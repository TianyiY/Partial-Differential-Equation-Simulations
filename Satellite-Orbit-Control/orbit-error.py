import math
import matplotlib.pyplot
import numpy

# These are used to keep track of the data we want to plot
h_array = []
error_array = []

total_time = 24. * 3600.  # s
g = 9.81  # m / s2
earth_mass = 5.97e24  # kg
gravitational_constant = 6.67e-11  # N m2 / kg2
radius = (gravitational_constant * earth_mass * total_time ** 2. / 4. / math.pi ** 2.) ** (1. / 3.)
speed = 2.0 * math.pi * radius / total_time


def acceleration(spaceship_position):
    vector_to_earth = - spaceship_position  # earth located at origin
    return gravitational_constant * earth_mass / numpy.linalg.norm(vector_to_earth) ** 3 * vector_to_earth


def calculate_error(num_steps):
    h=total_time/num_steps
    x=numpy.zeros([num_steps+1, 2])  # m
    v=numpy.zeros([num_steps+1, 2])  # m/s

    x[0,0]=radius
    v[0,1]=speed

    for step in range(num_steps):
        x[step+1]=x[step]+h*v[step]
        v[step+1]=v[step]+h*acceleration(x[step])

    error=numpy.linalg.norm(x[-1]-x[0])
    matplotlib.pyplot.scatter(h, error)
    # This is used for plotting
    h_array.append(h)
    error_array.append(error)
    return error


for num_steps in [200, 500, 1000, 2000, 5000, 10000]:
    error = calculate_error(num_steps)


def plot_me():
    axes = matplotlib.pyplot.gca()
    axes.set_xlabel('Step size in s')
    axes.set_ylabel('Error in m')
    matplotlib.pyplot.scatter(h_array, error_array)
    matplotlib.pyplot.xlim(xmin=0.)
    matplotlib.pyplot.ylim(ymin=0.)
    matplotlib.pyplot.show()


plot_me()
