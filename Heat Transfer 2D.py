import numpy as np
from matplotlib import pyplot

# Defined problem

diff = 110    #Heat Diffusivity
length = 50   #millimteres
time = 10      #econds
nodes = 20    #vector nodes


#Initialization
dx = length / nodes
dy = length / nodes

dt = min(dx**2 / (4*diff), dy**2 / (4*diff))
t_nodes = time//dt

vec = np.zeros((nodes, nodes)) + 20

vec[0, :] = 0
vec[-1, :] = 0

vec[:, 0] = 100
vec[:, -1] = 100

# Visualization
fig, axis = pyplot.subplots()

pcm = axis.pcolormesh(vec, cmap=pyplot.cm.jet, vmin=0, vmax=100)
pyplot.colorbar(pcm, ax=axis)

#Simulation
counter = 0

while counter < time:

    vecopy = vec.copy()

    for i in range(1, nodes-1):
        for j in range(1, nodes-1):

            dd_ux = (vecopy[i-1, j] - 2*vecopy[i, j] + vecopy[i+1, j])/dx**2
            dd_uy = (vecopy[i, j-1] - 2*vecopy[i, j] + vecopy[i, j+1])/dy**2

            vec[i, j] = dt * diff * (dd_ux + dd_uy) + vecopy[i, j]

    counter += dt
    pcm.set_array(vec)
    axis.set_title("Distribution at t: {:.3f} [s]\nAverage temperature: {:.2f} [c]"
                   .format(counter, np.average(vec)))
    pyplot.pause(0.01)



pyplot.show()



