import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import math

background_temperature = 273. # K
source_temperature = 1273. # K
Lambda = 25 # thermal conductivity (m^2 / s)
delta_x = 1 # m
grid_num = 100 # grids counts
delta_t = 1 # s
stop_time = 2000. # s
steps = int(stop_time / delta_t)
###################################################
heat_loss_time = 200.  # s
x_velocity = 0.15  # m / s
y_velocity = 0.12  # m / s
ignition_temperature = 561.  # K
flame_temperature = 1400.  # K
fuel_consumption_time = 1800.  # s
fuel_density=300.  # kg/m^2
heating_value = (flame_temperature - background_temperature) / ( heat_loss_time * fuel_density) * fuel_consumption_time  # K / (kg / m2)
# fuel distribution
fuel_slope = 0.4  # dimensionless
fuel_intercept_inner = 100.  # m
fuel_intercept_outer = 200.  # m
fuel_density_inner = fuel_density  # kg / m2
fuel_density_outer = fuel_density_inner*0.6  # kg / m2
length = 650.  # [-length, +length], meters
dx = 2. * length / grid_num
# Pick a time step below the threshold of instability
delta_t_fuel = 0.2 * dx ** 2 / Lambda  # s

def heat_conduction_explicit():
    # initial conditions
    temperatures_old = background_temperature * np.ones([grid_num, grid_num])
    for y in range(int(4.5 * grid_num / 10), int(5.5 * grid_num / 10)):
        for x in range(int(4.5 * grid_num / 10), int(5.5 * grid_num / 10)):
            temperatures_old[y, x] = source_temperature
    temperatures_new = np.copy(temperatures_old)
    for step in range(steps):
        for y in range(1, grid_num - 1):
            for x in range(1, grid_num - 1):
                temperatures=temperatures_old[y,x]
                temperatures_new[y, x]=temperatures+delta_t*Lambda/delta_x**2\
                                       *(temperatures_old[y, x-delta_x]+temperatures_old[y, x+delta_x]
                                         +temperatures_old[y-delta_x, x]+temperatures_old[y+delta_x, x]
                                         -4*temperatures)
        temperatures_old, temperatures_new = temperatures_new, temperatures_old
    return temperatures_old

def heat_conduction_implicit():
    temperatures_old = background_temperature * np.ones(grid_num)
    for i in range(int(4.5 * grid_num / 10), int(5.5 * grid_num / 10)):
        temperatures_old[i] = background_temperature
    temperatures_new = np.copy(temperatures_old)
    coeff = delta_t * Lambda / delta_x ** 2
    coefficients = np.zeros([grid_num, grid_num])
    for i in range(1, grid_num-1):
        coefficients[i, i]=1.+2.*coeff
    for i in range(0, grid_num-1):
        coefficients[i, i+1]=-coeff
        coefficients[i+1, i]=-coeff
    coefficients[0, 0]=1.+coeff
    coefficients[-1, -1]=1.+coeff
    for step in range(steps):
        temperatures_new=np.linalg.solve(coefficients, temperatures_old)
        temperatures_old, temperatures_new=temperatures_new, temperatures_old
    return temperatures_old


# Convert from integer grid positions to coordinates measured in meters
def grid2physical(i, j):
    return i * dx - length + 0.5 * dx, j * dx - length + 0.5 * dx

def heat_diffusion_with_HeatLoss_AirFlow_Combustion():
    temperatures_old = background_temperature * np.ones([grid_num, grid_num])  # K
    fuel_old = np.zeros([grid_num, grid_num])  # kg / m2
    for j in range(0, grid_num):
        for i in range(0, grid_num):
            x, y = grid2physical(i, j)
            temperatures_old[j][i]=(flame_temperature-background_temperature)*math.exp(-((x + 50.)**2+(y + 250.)**2)
                                    /(2. * 50. ** 2)) + background_temperature
            intercept=y-fuel_slope*x
            fuel=fuel_density_inner+(fuel_density_outer-fuel_density_inner)/(fuel_intercept_outer-fuel_density_inner)*\
                                    (intercept-fuel_intercept_inner)
            fuel_old[j][i]=min(fuel_density_inner, max(fuel_density_outer, fuel))
    temperatures_new = np.copy(temperatures_old)  # K
    fuel_new = np.copy(fuel_old)  # kg / m2
    steps = int(stop_time / delta_t_fuel)
    for step in range(steps):
        for j in range(1, grid_num - 1):
            for i in range(1, grid_num - 1):
                temperature = temperatures_old[j][i]
                fuel_density_change_rate=0.
                if temperature>=ignition_temperature:
                    fuel_density_change_rate=fuel_old[j][i]/fuel_consumption_time
                fuel_new[j][i]=fuel_old[j][i]-delta_t_fuel*fuel_density_change_rate
                temperatures_new[j][i]=temperature+delta_t_fuel*(Lambda/dx**2*(temperatures_old[j-1][i]+
                                       temperatures_old[j][i-1]+temperatures_old[j+1][i]+temperatures_old[j][i+1]
                                       -4.*temperature)+(background_temperature-temperature)/heat_loss_time
                                       -0.5/dx*(x_velocity*(temperatures_old[j+1][i]-temperatures_old[j-1][i]))
                                       +heating_value*fuel_density_change_rate)
        temperatures_old, temperatures_new = temperatures_new, temperatures_old
        fuel_old, fuel_new = fuel_new, fuel_old
    return temperatures_old, fuel_old

temperatures, fuel = heat_diffusion_with_HeatLoss_AirFlow_Combustion()

def plot_heat():
    axes = plt.gca()
    dimensions = [0, delta_x * grid_num, 0, delta_x * grid_num]
    plt.imshow(temperatures, cmap = matplotlib.cm.hot, origin = 'lower', extent = dimensions)
    plt.colorbar().set_label('Temperature')
    axes.set_xlabel('X Position')
    axes.set_ylabel('Y Position')
    plt.show()

def plot_fuel():
    dimensions = [-length, length, -length, length]
    axes = plt.gca()
    matplotlib.pyplot.imshow(fuel, cmap=matplotlib.cm.winter, origin='lower', extent=dimensions)
    matplotlib.pyplot.colorbar()
    axes.set_title('Density of Fuel')
    axes.set_xlabel('X Position')
    axes.set_ylabel('Y Position')
    plt.show()

plot_heat()
plot_fuel()