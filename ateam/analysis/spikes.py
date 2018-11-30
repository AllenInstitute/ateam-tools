import numpy as np 
# from bmtk.analyzer.spikes_loader import load_spikes
from bmtk.utils.spike_trains import SpikesFile
from ateam.sim.setup import SimManager

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