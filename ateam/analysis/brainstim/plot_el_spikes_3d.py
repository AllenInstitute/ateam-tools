######################################################
# Authors: Fahimeh Baftizadeh, Taylor Connington
# Date created: 9/1/2017
######################################################

import numpy as np
import matplotlib.pyplot as plt
import reobase_utils as ru
import analysis as ra
from mpl_toolkits.mplot3d import Axes3D
from .reobase_utils import StimType

# Get radius for soma sphere model
#from allensdk.core.cell_types_cache import CellTypesCache
#ctc = CellTypesCache(manifest_file='cell_types/manifest.json')
#morphology = ctc.get_reconstruction(cell_gid)
#soma = morphology.compartment_index[0]

#%% Plot

PLOT_SOMA = True

def plot(t,threshold=None):

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    if threshold is not None:
        yes = t[t['num_spikes'] >  threshold]
        no  = t[t['num_spikes'] <= threshold]
        ax.scatter(yes.x, yes.y, yes.z, c='g', s=50, edgecolors='none')
        ax.scatter(no.x, no.y, no.z, c='grey', s=50, edgecolors='none')
    else:
        line = ax.scatter(t.x, t.y, t.z, c=t.spikes.map(lambda x: x.size),
                          s=50, edgecolors='none')
        fig.canvas.mpl_connect('draw_event', lambda ev: line.changed())
        plt.colorbar(line)

    if PLOT_SOMA:
        u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
        r = 7  # soma['radius']
        x = r * np.cos(u) * np.sin(v)
        y = r * np.sin(u) * np.sin(v)
        z = r * np.cos(v)
        ax.plot_wireframe(x, y, z, color="r")

    plt.show()

#%%
    
def default():
    cell_gid = [313862022, 314900022, 320668879][1]
    t = ru.read_cell_tables(471077468, [str(x) for x in [30, 40]], stim_type=StimType.DC)
    t['theta'], t['phi'] = ra.spherical_coords(t)

    amp = -0.04
    sub = t[(t['amp'] == amp)] # only really makes sense if you filter by amp
    # sub = sub[(sub['phi'] > 3 * np.pi / 4) | (sub['phi'] < np.pi / 4)]

    plot(t, 0)

#default()