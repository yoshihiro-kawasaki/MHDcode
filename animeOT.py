import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation
import matplotlib.animation as animation
from matplotlib import gridspec
import re
import os
import linecache

dirname = "output/OrszagTang3"
nfile = 31

savename = "OT3.mp4"

fig = plt.figure(figsize=(12, 10))
ims = []

gridfile = dirname + "/grid"
x = np.array([float(i) for i in linecache.getline(gridfile, 1).split()])
y = np.array([float(i) for i in linecache.getline(gridfile, 2).split()])
nx = len(x)
ny = len(y)
X, Y = np.meshgrid(x, y)


def animate(idx):
    plt.gcf().clear()
    filename = dirname + "/st" + str(int(idx))
    is_file = os.path.isfile(filename)
    if (not is_file):
        return
    time = float(linecache.getline(filename, 1))
    print(idx, filename, time)
    linecache.clearcache()
    data = np.loadtxt(filename, dtype=float, skiprows=1).T
    rho = data[0]
    vx = data[1]
    vy = data[2]
    vz = data[3]
    p = data[4].reshape(nx, ny)
    bx = data[5]
    by = data[6]
    bz = data[7]
    img = plt.pcolormesh(X, Y, p)
    plt.title(" (t = {0:.2f})".format(time), fontsize=30)
    plt.xlim(x[0], x[-1])
    plt.tick_params(labelsize=30)
    #ims.append([img])

ani = animation.FuncAnimation(fig, animate, interval = 80, frames = nfile+1)
ani.save(savename, writer = "imagemagick")