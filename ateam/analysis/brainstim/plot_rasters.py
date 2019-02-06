######################################################
# Authors: Fahimeh Baftizadeh, Taylor Connington
# Date created: 9/1/2017
######################################################

import numpy as np
import matplotlib.pyplot as plt
import reobase_utils as r
from .reobase_utils import StimType


#%% find raster for a given el, amp

# by amp

def build_raster_data(t, slice_feature, value, y_feature):
    sub_df = t[t[slice_feature] == value]
    spike_times = np.array([]) # TODO preallocate
    ys = np.array([])
    dist = np.array([])

    for rid,row in sub_df.iterrows():
        feat = row[y_feature]
        num = row['num_spikes']
        spikes = row['spikes']
        d = row['distance']
        spike_times = np.concatenate((spike_times, spikes))
        ys = np.concatenate((ys, [feat]*num))
        dist = np.concatenate((dist, [d]*num))
        
    return spike_times, ys, dist

#%% electrode spike raster for specific amp

    
def plot_amp(t, amp):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    xs, ys, dist = build_raster_data(t, 'amp', amp, 'electrode')
    normed_dist = (dist - t.distance.min())/t.distance.max()
    ax.scatter(xs, ys, marker='o', s=10, c=plt.cm.viridis(normed_dist), edgecolors='none')
    
        
    ax.set_title('Spike raster for amp = {} (mA)'.format(amp))
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Electrode')
    plt.show()


def default():
    cell_gid = [313862022, 314900022, 320668879][1]
    stim_type = StimType.DC
    amp_range = [str(x) for x in range(10, 20, 10)]

    t = r.read_cell_tables(cell_gid, amp_range, stim_type=stim_type)
    for a in t.amp.unique():
        plot_amp(t, a)

print 'For an example try default() function'

# amp spike raster for specific electrode
#def plot_el(el):
#    xs, ys = build_raster_data(t, 'electrode', el, 'amp')
#    
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    ax.scatter(xs, ys,
#                      marker='o', s=10, edgecolors='none')
#
#    ax.set_title('Spike raster for el = {}'.format(el))
#    ax.set_xlabel('Time')
#    ax.set_ylabel('$I_{stim}$')
#    plt.show()

#plot_el(10)



