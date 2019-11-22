# Allen Institute Software License - This software license is the 2-clause BSD
# license plus a third clause that prohibits redistribution for commercial
# purposes without further permission.
#
# Copyright 2018. Allen Institute. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# 3. Redistributions for commercial purposes are not permitted without the
# Allen Institute's written permission.
# For purposes of this license, commercial purposes is the incorporation of the
# Allen Institute's software into anything for which you will charge fees or
# other compensation. Contact terms@alleninstitute.org for commercial licensing
# opportunities.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor
from allensdk.core.cell_types_cache import NwbDataSet
import numpy as np
import pandas as pd
import h5py
import scipy.signal
import warnings
from collections import namedtuple


DIDT_THRESH = 10 # Threshold for identifying input pulses by dI/dt
TRIPLE = 'Short Square - Triple'

def spikegroups(sweepex):
    # Group spikes by stimulus pulse
    index_onset = np.flatnonzero(np.diff(sweepex.i) > DIDT_THRESH)
    peak_index = sweepex.spike_feature("peak_index")
    spikes1 = np.flatnonzero((peak_index>index_onset[0]) & (peak_index<index_onset[1]))
    spikes2 = np.flatnonzero((peak_index>index_onset[1]) & (peak_index<index_onset[2]))
    spikes3 = np.flatnonzero(peak_index>index_onset[2] )
    return (spikes1, spikes2, spikes3)

def postspike_area(sweepex):
    v_slow = sweepex.spike_feature("slow_trough_v")
    if len(v_slow)<3:
        return np.nan

    vhalf_post = (v_slow[-1] + sweepex.spike_feature("fast_trough_v")[-1])/2
    term_index = sweepex.spike_feature("fast_trough_index")[-1]
    postspike_index = term_index + np.nonzero(sweepex.v[term_index:] < vhalf_post)[0][0]
    postspike_trace = sweepex.v[term_index: postspike_index]
    dt = sweepex.t[1] - sweepex.t[0]
    postspike_area = (np.sum(postspike_trace) - len(postspike_trace)*vhalf_post) *1000*dt
    return postspike_area

def postspike_v(sweepex):
    """Post-final-spike voltage"""
    return postspike_v_all(sweepex)[-1]

def postspike_v_all(sweepex, delta_t=0.01):
    """Post-spike voltage for all spikes"""
    spike_index = sweepex.spike_feature("peak_index")
    dt = sweepex.t[1] - sweepex.t[0]
    postspike_index = spike_index + int(delta_t/dt)
    post_v = try_index(sweepex.v, postspike_index)
    return post_v

def postspike_v_range(sweepex):
    spike_index = sweepex.spike_feature("peak_index")[-1]
    dt = sweepex.t[1] - sweepex.t[0]
    postspike_index_all = spike_index + (np.arange(0.008,0.019,0.001)/dt).astype(int)
    post_v_all = try_index(sweepex.v, postspike_index_all)
    return post_v_all

def postspike_v_relative(sweepex, delta_t=0.01, baseline_interval=0.01):
    post_v = postspike_v_all(sweepex, delta_t=0.01)[-1]
    index_onset = np.flatnonzero(np.diff(sweepex.i) > DIDT_THRESH)
    dt = sweepex.t[1] - sweepex.t[0]
    index_baseline = index_onset[0] - range(int(baseline_interval/dt))
    v_baseline = np.mean(sweepex.v[index_baseline])
    # initial_v = sweepex.v[index_onset[0]]
    return post_v - v_baseline

def vpost_delta(sweepex):
    post_v = postspike_v_all(sweepex, delta_t=0.005)
    i_first = spikegroups(sweepex)[0][0]
    return post_v[-1] - post_v[i_first]

def adp_delta(sweepex):
    list_vpost = sweepex.spike_feature("fast_trough_v")
    if len(list_vpost)<3:
        return np.nan

    list_vadp = sweepex.spike_feature("adp_v")
    list_adp = np.where(~np.isnan(list_vadp), list_vadp, list_vpost)
    return  list_adp[-1] - list_adp[0]

def amp(sweepex):
    index_onset = np.flatnonzero(np.diff(sweepex.i) > DIDT_THRESH)
    if index_onset.size>0:
        return sweepex.i[index_onset[0] + 1]
    else:
        return None

def frequency(sweepex):
    index_onset = np.flatnonzero(np.diff(sweepex.i) > DIDT_THRESH)
    if index_onset.size>0:
        t = sweepex.t
        freq = 1/(t[index_onset[1]] - t[index_onset[0]])
        return freq
    else:
        return None

def try_index(array, index):
    # Helper method for indexing that may go beyond end of array
    imax = index if np.isscalar(index) else max(index)
    if imax > len(array):
        array = np.pad(array, (0, imax-len(array)+1), mode='constant', constant_values=np.nan)
    return array[index]

def get_new_feature(sweepex, featurename, func):
    EphysSweepFeatureExtractor.process_new_sweep_feature(sweepex, featurename, func)
    return EphysSweepFeatureExtractor.sweep_feature(sweepex, featurename)

def sweepex_from_nwb_sweep(sweep):
    ii = sweep["index_range"]
    i = sweep["stimulus"][ii[0]: ii[1]+1] # in A
    v = sweep["response"][ii[0]: ii[1]+1]
    i *= 1e12 # to pA
    v *= 1e3 # to mV
    sampling_rate = sweep["sampling_rate"] # in Hz
    t = np.arange(0, len(v)) * (1.0 / sampling_rate)

    sweepex = EphysSweepFeatureExtractor(t=t, v=v, i=i)
    # sweepex.process_spikes()
    return sweepex

def sweepex_from_lims_nwb(nwb_path, sweepnum):
    # need to process nwb directly with h5py
    # NwbDataset class seems to interpret start/end times incorrectly for MET data 
    with h5py.File(nwb_path,'r') as nwb_file:
        return sweepex_from_lims_file(nwb_file, sweepnum)

def sweepex_from_lims_file(nwb_file, sweepnum):
    voltage_name = 'acquisition/timeseries/Sweep_%d/data' % sweepnum
    current_name = 'stimulus/presentation/Sweep_%d/data' % sweepnum
    rate_name = 'acquisition/timeseries/Sweep_%d/starting_time' % sweepnum
    # in mV
    voltage=nwb_file[voltage_name].value*1e3
    # in pA
    current=nwb_file[current_name].value*1e12
    # extracts the sampling rate for a given recording in Hz
    sampling_rate=nwb_file[rate_name].attrs['rate']
    # time is in sec

    tstart = 0
    idiff = np.diff(current)
    change_points = np.flatnonzero(idiff)

    num_sep = 1000 # padding after test pulse
    if len(change_points)>6: # more than 3 pulses indicates initial test pulse
        # trim trace to start after test pulse
        i_start = change_points[-7] + num_sep
        voltage = voltage[i_start:]
        current = current[i_start:]
    time = np.arange(0, len(voltage)) * (1.0 / int(sampling_rate))

    sweepex = EphysSweepFeatureExtractor(t=time, v=voltage, i=current, start=tstart)
    return sweepex

import logging
from functools import partial
from multiprocessing import Pool

import ateam.data.lims as lr
def get_tri_features_lims_single(cells_df, cell):
    limsreader = lr.LimsReader()
    nwb_path = cells_df.loc[cell].nwb_path
    sweeps_tri = limsreader.get_sweeps(cell, TRIPLE)
    del limsreader
    if not sweeps_tri:
        return None
    cell_data = dict()
    with h5py.File(nwb_path,'r') as nwb_file:
        try:
            for sweepnum in sweeps_tri:
                # check for null stim
                sweepex = sweepex_from_lims_file(nwb_file, sweepnum)
                tri_dict = get_tri_features_single(sweepex)
                if tri_dict:
                    cell_data[sweepnum] = tri_dict
        except Exception as detail:
            logging.warn("Exception when processing specimen {:d}".format(cell))
            logging.warn(detail)
    cell_df = pd.DataFrame.from_dict(cell_data, orient='index')
    return cell_df

def get_tri_features_lims(cells_df, run_parallel=True):
    cell_list = cells_df.index.values
    get_data_partial = partial(get_tri_features_lims_single, cells_df)
    if run_parallel:
        pool = Pool(processes=4)
        tri_data = pool.map(get_data_partial, cell_list)
    else:
        tri_data = map(get_data_partial, cell_list)
    # if stored as dict:
    # tri_df = pd.concat(tri_data, names=['cell', 'sweep'], sort=True)
    tri_df = pd.concat(tri_data, keys=cell_list, names=['cell', 'sweep'], sort=True)
    return tri_df

def has_triblip_lims(limsreader, cell_id):
    return len(limsreader.get_sweeps(cell_id, 'Short Square - Triple'))>0

def get_tri_features(ctc, cell_list):
    tri_data = dict()
    for cell in cell_list:
        cell_data = dict()

        dataset = ctc.get_ephys_data(cell)
        # sweepnum_list = [n for n in dataset.get_sweep_numbers() if
        #              dataset.get_sweep_metadata(n)['aibs_stimulus_name'] == TRIPLE]
        sweeps_df = pd.DataFrame(ctc.get_ephys_sweeps(cell))
        sweeps_tri = sweeps_df[sweeps_df.stimulus_name==TRIPLE]
        # Need to find the final amplitude to skip the threshold-finding sweeps
        amp_final = sweeps_tri.stimulus_absolute_amplitude.iloc[-1]
        sweeps_tri_final = sweeps_tri[ np.isclose(amp_final, sweeps_tri.stimulus_absolute_amplitude)]
        sweepnum_list = sweeps_tri_final.sweep_number

        for sweepnum in sweepnum_list:
            sweep = dataset.get_sweep(sweepnum)
            sweepex = sweepex_from_nwb_sweep(sweep)
            cell_data[sweepnum] = get_tri_features_single(sweepex)
        cell_df = pd.DataFrame.from_dict(cell_data, orient='index')
        tri_data[cell] = cell_df
    tri_df = pd.concat(tri_data, names=['cell', 'sweep'], sort=True)
    # tri_df = pd.concat(tri_data, keys=cell_list, names=['cell', 'sweep'], sort=True)

    tri_df.is_bursty = tri_df.is_bursty.astype('bool')
    tri_df.is_complete = tri_df.is_complete.astype('bool')
    tri_df["is_complete_noburst"] = tri_df.is_complete & ~tri_df.is_bursty
    return tri_df

def get_tri_features_single(sweepex):
    sweep_data = dict()
    sweep_data["amp"] = amp(sweepex)
    if sweep_data["amp"] is None:
        return None
    sweep_data["frequency"] = frequency(sweepex)
    try:
        sweepex.process_spikes()
    except ValueError:
        warnings.warn('Error processing spikes')
        return None

    spike_groups = spikegroups(sweepex)
    is_complete = all(spikes.size>0 for spikes in spike_groups)
    sweep_data["is_complete"] = is_complete
    is_bursty = any(spikes.size>1 for spikes in spike_groups)
    sweep_data["is_bursty"] = is_bursty
    if is_complete:
        sweep_data["postspike_v"] = postspike_v(sweepex)
        sweep_data["postspike_v_relative"] = postspike_v_relative(sweepex)
        # sweep_data["postspike_v_range"] = postspike_v_range(sweepex)
        # sweep_data["delta_adp"] = adp_delta(sweepex)
        # sweep_data["delta_vpost"] = vpost_delta(sweepex)
    return sweep_data

def vpost_trend(x):
    freq_upper = 100
    freq_lower = 65
    return x[x.frequency>freq_upper].postspike_v.mean() \
        - x[(10<x.frequency)&(x.frequency<freq_lower)].postspike_v.mean()

def process_tri_trend(tri_df):
    groups_cell = tri_df.groupby('cell')
    series_trend = groups_cell.apply(vpost_trend).dropna()
    series_trend.name = "trend_vpost"
    return series_trend

# Sigmoid functions
import scipy.optimize as opt
def df_mask_outliers(df, factor=4):
    lq = df.quantile(0.25)
    uq = df.quantile(0.75)
    iqr = uq-lq
    df = (df > uq + factor*iqr) | (df < lq - factor*iqr)
    return df

def fsigmoid(f, v0, dv, f0, k):
    return v0 + dv / (1.0 + np.exp(-k*(f-f0)))

def fit_sigmoid(celldat, feature='postspike_v'):
    min_sweeps = 10
    f0_bound = [50, 150]
    k_bound = [0.01, 0.5]
    dv_bound = [-10, 10]
    celldat = celldat[celldat.is_complete & ~celldat.outlier]
    if len(celldat.is_complete) < min_sweeps:
        raise Exception('Not enough triblip sweeps')
    xdata = np.array(celldat.frequency)
    ydata = np.array(celldat[feature])
    v_bound = [min(ydata), max(ydata)]
    bounds = list(zip(v_bound, dv_bound, f0_bound, k_bound))
    p0 = [np.mean(ydata), 0, 80, 0.05]
    
    popt, pcov = opt.curve_fit(fsigmoid, xdata, ydata, method='trf', bounds=bounds, p0=p0)
    triblip_dv = fsigmoid(150, *popt) - fsigmoid(50, *popt)
    return list(popt)+[triblip_dv]
            
def process_tri_sigmoid(tri_df, relative=False):
    feature = 'postspike_v' if not relative else 'postspike_v_relative'
    tri_df['outlier'] = tri_df.groupby('cell')[feature].transform(df_mask_outliers)
    cells = tri_df.index.unique(level='cell')
    fit_data = {}
    for cell in cells:
        try:
            celldat = tri_df.loc[cell]
            popt = fit_sigmoid(celldat, feature=feature)
            fit_data[cell] = popt
        except Exception as detail:
            warnings.warn("Exception when processing specimen {:d}".format(cell))
            warnings.warn(detail)
    fit_df = pd.DataFrame.from_dict(fit_data, orient='index', columns=['v0','dv','f0','k','triblip_dv'])
    return fit_df

# Population analysis
def process_groups(combined_df, feature_col, group_col):
    popdata = combined_df[feature_col]
    popstats = pop_stats_robust(popdata)
    def frac_str(col):
        nsig = (col>popstats.cutoff).sum()
        ntot = col.count()
        return '%d/%d' % (nsig, ntot)
    def zscore(col):
        return (col.mean() - popstats.mean)/np.sqrt(popstats.std/col.count())
    def mw_test(col):
        import scipy.stats as stats
        u, p = stats.mannwhitneyu(popdata, col, alternative='less')
        sign = np.sign(u - len(popdata)*len(col)/2)
        return p
    return combined_df.groupby(group_col)[feature_col].agg([frac_str, zscore, mw_test])

def pop_stats_robust(series_trend, cutoff_z=3):
    data = series_trend.values
    robust_stats = namedtuple("robust_stats", ["mean", "std", "cutoff"])
    trend_med = np.median(data)
    trend_trimean = (np.percentile(data, 25) + 2*trend_med + np.percentile(data, 75))/4
    trend_std = 1.4826 * np.median(np.abs(data - trend_med))
    cutoff = trend_trimean + cutoff_z*trend_std
    return robust_stats(trend_trimean, trend_std, cutoff)

# TODO: remove
def load_ttype_df(shiny_path):
    ttype_df = pd.read_feather(shiny_path)
    ttype_df.spec_id_label.replace("ZZ_Missing", None, inplace=True)
    ttype_df['spec_id_label'] = pd.to_numeric(ttype_df['spec_id_label'])
    ttype_df.set_index('spec_id_label', inplace=True)
    return ttype_df

def process_features_for_shiny(tri_df, shiny_path, feather_path, csv_path):
    series_trend = process_tri_trend(tri_df)
    ttype_df = load_ttype_df(shiny_path)
    combined_df = ttype_df.join(series_trend, how='inner')
    popstats = pop_stats_robust(series_trend)
    # frac = lambda col: frac_str(col, popstats)
    # z = lambda col: zscore(col, popstats)
    # p = lambda col: mw_test(col, series_trend)
    # results = combined_df.groupby('cluster_label').trend_vpost.agg([frac, z, p])
    results = process_groups(combined_df, 'cluster_label')
    results['bar_col'] = results.zscore
    results.to_csv(csv_path)

    feather_df = pd.read_feather(shiny_path)
    feather_df.columns = feather_df.columns.astype(str)
    feather_df.spec_id_label.replace("ZZ_Missing", None, inplace=True)
    feather_df.spec_id_label = pd.to_numeric(feather_df.spec_id_label)
    tri_df_response = series_trend > popstats.cutoff
    tri_df_response.name = 'ephys_triple_pulse_responsive_label'
    feather_df = feather_df.join(tri_df_response)
    # ttype_df.reset_index(inplace=True)
    feather_df.to_feather(feather_path)

    return combined_df