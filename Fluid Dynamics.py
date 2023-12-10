import numpy as np
from matplotlib import pyplot

PLOT_EVERY = 10 # skip some to make this run smoothly

def distance(x1, y1, x2, y2):
    return np.sqrt((x2-x1)**2 + (y2-y1)**2)

def main():
    Nx = 400 # x-axis
    Ny = 100 # y-axis
    tau = 0.53 #timescale
    Nt = 4000 #iterations

    # lattice speed and weights
    Nl = 9 #
    cxs = np.array([0, 0, 1, 1,  1,  0, -1, -1, -1])
    cys = np.array([0, 1, 1, 0, -1, -1, -1,  0,  1])
    weights = np.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36])

    F = np.ones((Ny, Nx, Nl)) + 0.01 * np.random.randn(Ny, Nx, Nl)
    F[:, :, 3] = 2.3 # for every y and x, 3rd node. This is to start moving rightward

    cylinder = np.full((Ny, Nx), False)

    for y in range(0, Ny):
        for x in range(0, Nx):
            if (distance(Nx//4, Ny//2, x, y) < 13):
                cylinder[y][x] = True
    
    # main loop
    for iteration in range(Nt):
        print(iteration)

        F[:, -1, [6, 7, 8]] = F[:, -2, [6, 7, 8]]
        F[:,  0, [2, 3, 4]] = F[:,  1, [2, 3, 4]]

        #all velocities to neighboring node lattices
        for i, cx, cy in zip(range(Nl), cxs, cys):
            F[:, :, i] = np.roll(F[:, :, i], cx, axis = 1)
            F[:, :, i] = np.roll(F[:, :, i], cy, axis = 0)

        # invert velocity in boundries
        boundryF = F[cylinder, :]
        boundryF = boundryF[:, [0, 5, 6, 7, 8, 1, 2, 3, 4]] #setting opposite velocity on collision

        # fluid variables
        rho = np.sum(F, 2) # summation for density
        ux = np.sum(F * cxs, 2) / rho
        uy = np.sum(F * cys, 2) / rho

        # set velocities inside cylinder to 0
        F[cylinder, :] = boundryF
        ux[cylinder] = 0
        uy[cylinder] = 0

        # collisions
        Feq = np.zeros(F.shape)
        for i, cx, cy, w in zip(range(Nl), cxs, cys, weights):
            # BGK approximation
            Feq[:, :, i] = rho * w * (
                1 + 3 * (cx*ux + cy*uy) + 9 * (cx*ux + cy*uy)**2 /2 - 3 *(ux**2 + uy**2) / 2
            )

        F = F + -(1/tau) * (F-Feq)

        if(iteration % PLOT_EVERY == 0):
            pyplot.imshow(np.sqrt(ux**2 + uy**2))
            pyplot.pause(0.01)
            pyplot.cla()

main()