import numpy as np
from matplotlib import pyplot

# Defined problem

diff = 110    #Heat Diffusivity
length = 50   #millimteres
time = 10      #econds
nodes = 20    #vector nodes


#Initialization
dx = length / nodes
dt = 0.5 * dx**2 / diff
t_nodes = time//dt

vec = np.zeros(nodes) + 20

vec[0] = 100
vec[-1] = 100

# Visualization
fig, axis = pyplot.subplots()

pcm = axis.pcolormesh([vec], cmap=pyplot.cm.jet, vmin=0, vmax=100)
pyplot.colorbar(pcm, ax=axis)
axis.set_ylim([-2,3])

#Simulation
counter = 0

while counter < time:

    vecopy = vec.copy()

    for i in range(1, nodes-1):

        vec[i] = dt * diff * (vecopy[i-1] - 2 * vec[i] + vecopy[i+1]) / dx**2 + vec[i]

    counter += dt
    pcm.set_array([vec])
    axis.set_title("Distribution at t: {:.3f} [s]\nAverage temperature: {:.2f} [c]"
                   .format(counter, np.average(vec)))
    pyplot.pause(0.01)



pyplot.show()



