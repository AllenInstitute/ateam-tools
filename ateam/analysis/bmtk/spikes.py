import numpy as np 
import h5py
import warnings
# from bmtk.analyzer.spikes_loader import load_spikes
from bmtk.utils.spike_trains import SpikesFile
from .nodes import create_node_table
from ateam.sim.setup import SimManager
import bmtk.analyzer.visualization.spikes as vs
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec

def get_spike_file(config_file):
    sm = SimManager(config_file)
    return sm.spikes_file

def get_rates_config(config_file):
    sm = SimManager(config_file)
    spikes = SpikesFile(sm.spikes_file)
    # gids, counts = np.unique(spikes_dict['gid'], return_counts=True)
    gids = spikes.gids
    counts = np.array([len(spikes.get_spikes(gid)) for gid in gids])
    rates = counts / sm.sim_time * 1000.
    return dict(zip(gids, rates))


def plot_spikes_rates(config_file, netname, group_key=None, exclude=[], color_dict=None, cmap='hsv', bins=100, fig=None, tlim=None):
    sm = SimManager(config_file)
    tlim = tlim or [0, sm.sim_time]
    fig = fig or plt.figure()
    nodes_df = create_node_table(sm.nodes_file(netname), sm.node_types_file(netname), group_key=group_key, exclude=exclude)

    # TODO: Uses utils.SpikesReader to open
    spikes_h5 = h5py.File(sm.spikes_file, 'r')
    spike_gids = np.array(spikes_h5['/spikes/gids'], dtype=np.int)
    spike_times = np.array(spikes_h5['/spikes/timestamps'], dtype=np.float)
    if len(spike_times)==0:
        warnings.warn("No spikes to plot")
        return

    if group_key is not None:
        # TODO: order by gid if groups are separate??
        if group_key not in nodes_df:
            raise Exception('Could not find column {}'.format(group_key))
        groupings = nodes_df.groupby(group_key)

        if color_dict is not None:
            groups_sub = []
            color_map = []
            for key, group in groupings:
                color = color_dict.get(key)
                if color:
                    groups_sub.append((key, group))
                    color_map.append(color)
        else:
            n_colors = nodes_df[group_key].nunique()
            color_norm = colors.Normalize(vmin=0, vmax=(n_colors-1))
            scalar_map = cmx.ScalarMappable(norm=color_norm, cmap=cmap)
            color_map = [scalar_map.to_rgba(i) for i in range(0, n_colors)]
    else:
        groupings = [(None, nodes_df)]
        color_map = ['blue']

    #marker = '.' if len(nodes_df) > 1000 else 'o'
    marker = 'o'

    # Create plot
    ngroups = len(color_map)
    gs = gridspec.GridSpec(ngroups+1, 1, height_ratios=[7] + ngroups*[1])

    ax1 = plt.subplot(gs[0])
    gid_min = nodes_df.index.min()
    gid_max = nodes_df.index.max()
    for color, (group_name, group_df) in zip(color_map, groupings):
        gids_group = group_df.index
        indexes = np.in1d(spike_gids, gids_group)
        ax1.scatter(spike_times[indexes], spike_gids[indexes], marker=marker, facecolors=color, label=group_name, lw=0, s=5)

    #ax1.set_xlabel('time (s)')
    ax1.axes.get_xaxis().set_visible(False)
    ax1.set_ylabel('cell_id')
    ax1.set_xlim(tlim)
    ax1.set_ylim([gid_min, gid_max])
    plt.legend(markerscale=2, scatterpoints=1)

    for i, group in enumerate(groupings):
        (group_name, group_df) = group
        color = color_map[i]
        gids_group = group_df.index
        indexes = np.in1d(spike_gids, gids_group)

        ax = plt.subplot(gs[i+1])
        plot_pop_rate(ax, spike_times[indexes], spike_gids[indexes], color=color, bins=bins, tlim=tlim)
        ax.axes.get_xaxis().set_visible(False)

    fig.text(0.05, 0.3, 'Firing rate (Hz)', rotation='vertical')

    ax.axes.get_xaxis().set_visible(True)
    ax.set_xlabel('time (ms)')

def plot_pop_rate(ax, spike_times, spike_gids, bins=100, tlim=None, **kwargs):
    n = len(np.unique(spike_gids))
    hist, bin_edges = np.histogram(spike_times, bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
    rates = hist / np.diff(bin_edges) * 1000 / n
    ax.plot(bin_centers, rates, **kwargs)
    # ax.set_xlabel('time (ms)')
    if tlim:
        ax.set_xlim(tlim)
    ax.set_ylim(bottom=0)

def plot_spikes_rates_traces(config_file, netname, gids=None, group_key=None, exclude=[], color_dict=None, cmap='hsv', bins=100, fig=None):
    from .cell_vars import plot_v
    sm = SimManager(config_file)
    fig = fig or plt.figure()
    nodes_df = create_node_table(sm.nodes_file(netname), sm.node_types_file(netname), group_key=group_key, exclude=exclude)

    # TODO: Uses utils.SpikesReader to open
    spikes_h5 = h5py.File(sm.spikes_file, 'r')
    spike_gids = np.array(spikes_h5['/spikes/gids'], dtype=np.int)
    spike_times = np.array(spikes_h5['/spikes/timestamps'], dtype=np.float)
    if len(spike_times)==0:
        warnings.warn("No spikes to plot")
        return

    if group_key is not None:
        # TODO: order by gid if groups are separate??
        if group_key not in nodes_df:
            raise Exception('Could not find column {}'.format(group_key))
        groupings = nodes_df.groupby(group_key)

        if color_dict is not None:
            groups_sub = []
            color_map = []
            for key, group in groupings:
                color = color_dict.get(key)
                if color:
                    groups_sub.append((key, group))
                    color_map.append(color)
        else:
            n_colors = nodes_df[group_key].nunique()
            color_norm = colors.Normalize(vmin=0, vmax=(n_colors-1))
            scalar_map = cmx.ScalarMappable(norm=color_norm, cmap=cmap)
            color_map = [scalar_map.to_rgba(i) for i in range(0, n_colors)]
    else:
        groupings = [(None, nodes_df)]
        color_map = ['blue']

    #marker = '.' if len(nodes_df) > 1000 else 'o'
    marker = 'o'

    # Create plot
    ngroups = len(color_map)
    gs = gridspec.GridSpec(ngroups+2, 1, height_ratios=[2] + [7] + ngroups*[1])
    imain = 1

    ax = plt.subplot(gs[0])
    plot_v( config_file, gids=gids, colors=color_map)
    ax.axis('off') 
    ax.axis(ymin=-80, ymax=50, xmin=0, xmax=sm.sim_time)

    ax1 = plt.subplot(gs[imain])
    gid_min = nodes_df.index.min()
    gid_max = nodes_df.index.max()
    for color, (group_name, group_df) in zip(color_map, groupings):
        gids_group = group_df.index
        indexes = np.in1d(spike_gids, gids_group)
        ax1.scatter(spike_times[indexes], spike_gids[indexes], marker=marker, facecolors=color, label=group_name, lw=0, s=5)

    #ax1.set_xlabel('time (s)')
    ax1.axes.get_xaxis().set_visible(False)
    ax1.set_ylabel('cell_id')
    ax1.set_xlim([0, max(spike_times)])
    ax1.set_ylim([gid_min, gid_max])
    plt.legend(markerscale=2, scatterpoints=1)

    for i, group in enumerate(groupings):
        (group_name, group_df) = group
        color = color_map[i]
        gids_group = group_df.index
        indexes = np.in1d(spike_gids, gids_group)

        ax = plt.subplot(gs[imain+i+1])
        plot_pop_rate(ax, spike_times[indexes], spike_gids[indexes], color=color, bins=bins)
        ax.axes.get_xaxis().set_visible(False)

    fig.text(0.05, 0.3, 'Firing rate (Hz)', rotation='vertical')

    ax.axes.get_xaxis().set_visible(True)
    ax.set_xlabel('time (ms)')