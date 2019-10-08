"""Analysis functions that operate directly on instantiated cells, without running simulations.
"""
from neuron import h
import bmtk.simulator.bionet.nrn as nrn
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from collections import defaultdict
import scipy.stats
import warnings

from ateam.sim.setup.sim_manager import create_singlecell_default
from ateam.sim.run.run_bionet import create_sim
try:
    from ateam.data.lims import node_dict_from_lims                     
except ImportError:
    pass

def sim_impedances(cell_id, sim_folder, freqs=[0,100], **sim_args):
    sim = sim_singlecell_transient(cell_id, sim_folder, **sim_args)
    imp = [morph_impedance(sim, freq=freq, detailed=True) for freq in freqs]
    morph = morph_props(sim)
    sim = None
    nrn.reset()
    return imp, morph

def sim_singlecell_transient(cell_id, sim_folder, active=True, sim_time=500, from_lims=False, directed=True):
    template = "/allen/aibs/mat/tmchartrand/bmtk_networks/biophys_components_shared/fast_config.json"
    node_dict = {}
    if from_lims:
        try:
            node_dict = node_dict_from_lims(cell_id, model_type='peri')
        except (NameError, UserWarning):
            warnings.warn("Failed to load cell from LIMS.")
            raise
    sm = create_singlecell_default(cell_id, sim_folder, sim_time=sim_time, active=active, directed=directed, config_template=template, node_dict=node_dict)
    sim = create_sim(sm.config_path)
    sim.run()
    return sim

def morph_props(sim, gid=0):
    cell = sim.net.get_cell_gid(gid)
    seg_props = cell.morphology.seg_prop
    nsegs = len(seg_props['x'])
    seg_coords = cell._seg_coords
    for name in ['p0', 'p1']:
        seg_props.update({name+'_x': seg_coords[name][0], 
                          name+'_y': seg_coords[name][1], 
                          name+'_z': seg_coords[name][2]})
    seg_props.update(seg_index=range(nsegs), gid=nsegs*[gid])
    return seg_props

def morph_impedance(sim, freq=0, gid=0, detailed=True, all_segs=True):
    cell = sim.net.get_cell_gid(gid)
    soma = cell.hobj.soma[0]
    imp_data = defaultdict(list)
    x_soma = 0.5
    imp = h.Impedance()
    imp.loc(x_soma, sec=soma)
    imp.compute(freq, detailed, sec=soma)
    rin_soma = imp.input(x_soma, sec=soma)
    # phase_soma = imp.input_phase(x_soma, sec=soma)
    if all_segs:
        for seg in cell.get_segments():
            _save_impedance(imp_data, imp, rin_soma, seg.sec, x=seg.x)
    else:
        for sec in cell.get_sections():
            _save_impedance(imp_data, rin_soma, imp, sec, x=0.5)
    for key, val in imp_data.items():
        imp_data[key] = np.array(val)
    imp_data['freq'] = freq
    return imp_data

def _save_impedance(imp_data, imp, rin_soma, sec, x):
    imp_data['rin'].append(imp.input(x, sec=sec))
    imp_data['transfer'].append(imp.transfer(x, sec=sec))
    imp_data['transfer_phase'].append(imp.transfer_phase(x, sec=sec))
    imp_data['ratio_out'].append(imp.transfer(x, sec=sec)/rin_soma)
    imp_data['ratio'].append(imp.ratio(x, sec=sec))

def points_to_img(values, xlims=None, ylims=None, n=100):
#     values shape: (# of dims, # of data)
    if xlims is None:
        xlims = [np.min(values[0,:]), np.max(values[0,:])]
    if ylims is None:
        ylims = [np.min(values[1,:]), np.max(values[1,:])]
    X, Y = np.mgrid[xlims[0]:xlims[1]:n*1j, ylims[0]:ylims[1]:n*1j]

    kernel = scipy.stats.gaussian_kde(values, bw_method=None)
    positions = np.vstack([X.ravel(), Y.ravel()])
    Z = np.reshape(kernel(positions).T, X.shape)
    return Z

def morph_impedance_coords(imp, morph, imp_prop='ratio_out', morph_prop='dist'):
    x = morph[morph_prop]
    y = imp[imp_prop]
    return x, y

def morph_impedance_img(imp, morph, xlims=[0,1500], **coords_args):
    x, y = morph_impedance_coords(imp, morph, **coords_args)
    Z = points_to_img(np.array([x,y]), xlims=xlims, ylims=[0,1], n=100)
    return Z

def morph_impedance_features(imp, morph, **coords_args):
    x, y = morph_impedance_coords(imp, morph, **coords_args)
    features = {
        'xmean': np.mean(x),
        'xvar': np.var(x),
        'xskew': scipy.stats.skew(x),
        'xkurtosis': scipy.stats.kurtosis(x),
        'ymean': np.mean(y),
        'yvar': np.var(y),
        'yskew': scipy.stats.skew(y),
        'ykurtosis': scipy.stats.kurtosis(y),
        'crosscorr': scipy.stats.pearsonr(x,y)[0]
    }
    return features


def get_impedance_img_single(cell_id, sim_folder, freq=100, sim_time=200):
    sim = sim_singlecell_transient(cell_id, sim_folder, active=True, sim_time=sim_time)
    imp = morph_impedance(sim, freq=freq, detailed=True)
    morph_dict = morph_props(sim)
    x = morph_dict['dist']
    y = imp['ratio_out']
    Z = points_to_img(np.array([x,y]), xlims=[0,1500], ylims=[0,1], n=100)
    return Z

