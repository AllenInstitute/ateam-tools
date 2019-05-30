
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import pandas as pd
import numpy as np
# In order: soma, axon, basal, apical
# numpy array to allow indexing with list
sec_colors = np.array(['', 'k', 'tab:blue', 'tab:red', 'tab:orange'])

# TODO: make compatible with AllenSDK Morphology class to use for both swc and BMTK morph

def plot_prop_dist_scatter(prop, seg_props, **kwargs):
    plot = plt.scatter
    dist = seg_props['dist']
    segtype = seg_props['type']
    # plot(prop, dist, color=sec_colors[segtype])
    plot(dist[segtype==2], prop[segtype==2], color='tab:blue', label='axon', **kwargs)
    plot(dist[segtype==3], prop[segtype==3], color='tab:red', label='dendrites', **kwargs)
    plot(dist[segtype==4], prop[segtype==4], color='tab:orange', label='apical dendrites', **kwargs)
    plot(dist[segtype==1], prop[segtype==1], color='k', label='soma', s=50, **kwargs)
    plt.xlabel("Distance from soma (um)")
    plt.legend()

def plot_segments_cmap(segments, color, cmap=None, centered=False, ax=None, **kwargs):
    ax = ax or plt.gca()
    plt.sca(ax)
    norm = None
    cmap = cmap or ('seismic' if centered else 'viridis')
    if centered:
        val = np.max(np.abs(color))
        norm = plt.Normalize(-val, val)
    lc = LineCollection(segments, array=np.array(color), cmap=cmap, norm=norm, **kwargs)
    line = ax.add_collection(lc)
    cbar = plt.colorbar(line, ax=ax)
    plt.autoscale()
    plt.xticks([])
    plt.yticks([])
    plt.axis('equal')

def plot_prop_morph_xy(prop, morph, cmap=None, centered=False, ax=None):
    segments = np.array([[morph['p0_x'],morph['p1_x']], [morph['p0_y'],morph['p1_y']]]).T
    plot_segments_cmap(segments, prop, cmap=cmap, centered=centered, ax=ax)

def plot_morph_yina(df):
    soma_idx = df['type']==1
    axon_idx = df['type']==2
    apical_idx = df['type']==4
    basal_idx= df['type']==3
    fig, axes = plt.subplots(1,2, figsize=(6,3))
    axes[0].plot([df['p0_x'],df['p1_x']],[df['p0_y'],df['p1_y']],color='k',zorder=-10)
    axes[0].set_ylabel('y')
    axes[0].set_xlabel('x')
    axes[0].plot([df['p0_x'][soma_idx],df['p1_x'][soma_idx]],[df['p0_y'][soma_idx],df['p1_y'][soma_idx]],color='darkred',linewidth='4')
    axes[0].plot([df['p0_x'][axon_idx],df['p1_x'][axon_idx]],[df['p0_y'][axon_idx],df['p1_y'][axon_idx]],color='cyan',zorder=10)
    simpleaxis(axes[0])

    axes[1].plot([df['p0_z'],df['p1_z']],[df['p0_y'],df['p1_y']],color='k')  
    axes[1].set_ylabel('y')
    axes[1].set_xlabel('z')    
    axes[1].plot([df['p0_z'][soma_idx],df['p1_z'][soma_idx]],[df['p0_y'][soma_idx],df['p1_y'][soma_idx]],color='darkred',linewidth='4')
    axes[1].plot([df['p0_z'][axon_idx],df['p1_z'][axon_idx]],[df['p0_y'][axon_idx],df['p1_y'][axon_idx]],color='cyan',zorder=10)
    simpleaxis(axes[1])
    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.4,wspace=0.3)

def simpleaxis(ax):
    #Hide the right and top spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
