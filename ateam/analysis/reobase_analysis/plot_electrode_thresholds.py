######################################################
# Authors: Fahimeh Baftizadeh, Taylor Connington
# Date created: 9/1/2017
######################################################

import matplotlib.pyplot as plt
import reobase_utils as r
import analysis as ra
import numpy as np
import reobase_analysis.tchelpers as tc

cell_gid = [313862022, 314900022, 320668879][2]
t = r.read_cell_tables(cell_gid, [str(x) for x in range(10,60,10)], stim_type='dc_lgn_poisson')
t['theta'], t['phi'] = ra.spherical_coords(t)

#%% Find threshold
## Only really makes sense for data with absent or subthreshold external input
thresholds = ra.find_thresholds(t)

el_th = []
for el,thresh in thresholds.items():
    if thresh is not None:
        el_th.append((el, thresh[1]))

[els, thresh] = zip(*el_th)

#%%
def get_phi(el):
    return t[t['electrode'] == el].phi[0]

def get_theta(el):
    return t[t['electrode'] == el].theta[0]


def plot_threshold_angle(els, thresh):

    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)
    line = ax.scatter([get_phi(el) for el in els], [get_theta(el) for el in els],
                c=thresh, s=14)
    plt.colorbar(line)
    ax.set_title('threshold at location in spherical coordinates')
    ax.set_xlabel('theta')
    tc.set_phi_ticks(ax)
    tc.set_theta_ticks(ax, 'y')
    plt.show()

#%%
def plot_phi():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter([get_phi(el) for el in els], thresh, s=4)
    ax.set_title('threshold vs angle $\phi$ to z-axis')
    tc.set_phi_ticks(ax)
    plt.show()

def plot_theta():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter([get_theta(el) for el in els], thresh, s=4)
    ax.set_title('threshold vs angle $\Theta$ to x-axis')
    tc.set_theta_ticks(ax)
    plt.show()

#%%
def table_covariance(t):
    cov = t[['x','y','z','phi','theta','distance','amp','num_spikes']].cov().as_matrix()
    np.fill_diagonal(cov, 0)
    
    line = plt.matshow(cov)
    plt.colorbar(line)
    plt.show()

#table_covariance(t)

