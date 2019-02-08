######################################################
# Authors: Fahimeh Baftizadeh, Taylor Connington
# Date created: 9/1/2017
######################################################

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import numpy as np
import matplotlib
import pandas as pd


def plot_3d_colorbar(X, Y, Z, Z1, cbar_min, cbar_max, figsize=(15,15), ax=None):

    fig = plt.figure(figsize=figsize)
    if not ax:
        ax = fig.gca(projection='3d')

    ax.tick_params(axis='x', which='major', pad=0)
    ax.tick_params(axis='z', which='major', pad=0)

    norm = matplotlib.colors.Normalize(vmin=cbar_min,vmax = cbar_max)

    surf1 = ax.plot_surface(X, Y, Z, rstride=1, cstride=1,facecolors=plt.cm.jet(norm(Z)),
                       linewidth=1, antialiased=True)
    surf2 = ax.plot_surface(X, Y, Z1, rstride=1, cstride=1, facecolors=plt.cm.jet(norm(Z1)),
                       linewidth=1, antialiased=True)

    m = cm.ScalarMappable(cmap=plt.cm.jet, norm=norm)
    m.set_array([])
    #cbar=plt.colorbar(m)
    #cbar.ax.tick_params(labelsize=20)
    ax.tick_params(labelsize=35)
    plt.gca().invert_yaxis()
    # ax.set_xticks([1,8,30,60,100])
    # ax.set_xticklabels([""])
    return ax
