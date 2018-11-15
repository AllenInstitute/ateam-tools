######################################################
# Authors: Fahimeh Baftizadeh
# Date created: 4/1/2018
######################################################

import glob
import numpy as np
import pandas as pd
import itertools
from scipy.signal import hilbert, chirp
import bmtk.analyzer.network_utils as nu
from  bmtk.analyzer.network_utils import StimType, ModelType, InputType
from bmtk.simulator.bionet.modules.xstim_waveforms import stimx_waveform_factory
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor

# import reobase_analysis.reobase_utils as ru
# import reobase_analysis.sin_utils as su
# from isee_engine.bionet.stimxwaveform import stimx_waveform_factory
# from isee_engine.bionet.stimxwaveform import iclamp_waveform_factory
# #import cProfile, pstats, io
# #pr = cProfile.Profile()
# #pr.enable()

"""
Script to build h5 files containing organized run info for all runs of a particular amplitude
Will include all electrodes. Has no overwrite protection.

Will attempt to include vm information if the stim type is DC. This is only done for DC because it only makes sense for
a constant input without external input. This allows us to look at the subthreshold response.
If the cell is active (> 1 spikes) then the post-stim vm value is NaN
"""

def build_table(network_id, input_type, model_type, stim_type, input_amp,  trial, include_spike_phase = True, saved_data= False, input_freq=None):

    include_sin = stim_type == "sin"
    out_dir = nu.get_output_dir(input_type = input_type, stim_type=stim_type, model_type= model_type, saved_data = saved_data, network_id=network_id)
    run_dir = nu.get_base_dir(saved_data)
    additional_cols = nu.resolve_additional_cols(include_sin, include_spike_phase)
    amp = input_amp
    freq = input_freq

    print 'Build data table for ex_amp = {}...'.format(amp)
    output_filename = nu.get_output_filename(amp, trial, freq)
    output = nu.concat_path(out_dir, output_filename)
    table = nu.build_df(additional_cols)
    amp = int([x for x in output_filename.split('_') if x.startswith('amp')][0][3:]) * 0.000001
    fq = int([x for x in output_filename.split('_') if x.startswith('freq')][0][4:]) if "freq" in output_filename else None
    config_path = nu.get_config_resolved_path(stim_type, output, amp, fq)
    cvh5 = nu.get_cv_files(output)
    nodesh5 = nu.get_nodes_files(network_id, run_dir)
    spikesh5 = nu.get_spikes_file(output)
    conf = nu.get_json_from_file(config_path)
    electrode_position_file = nu.get_dir_root(saved_data) + '/' + '/'.join(conf['inputs']['Extracellular_Stim']['positions_file'].split('/')[3:])
    el_xyz = nu.get_electrode_xyz(electrode_position_file).flatten()
    gids = cvh5['mapping/gids'].value
    node_ids = nu.get_node_id(nodesh5, network_id)
    node_type_ids = nu.get_node_type_id(nodesh5, network_id)
    node_group_ids = nu.get_node_group_id(nodesh5, network_id)
    node_group_index = nu.get_node_group_index(nodesh5, network_id)

    print "Wrting the basic table"
    for gid in gids:
        # vm_rest = np.NaN
        node_type_id = node_type_ids[gid]
        cell_group_id = node_group_ids[gid]
        cell_group_index = node_group_index[gid]
        cell_xyz = nu.get_cell_xyz(nodesh5, network_id, cell_group_id)[cell_group_index]
        el_dist = np.linalg.norm(el_xyz - cell_xyz)
        run_id = nu.resolve_run_id(gid, amp, fq)
        spikes = nu.get_cell_spikes(spikesh5, gid)

        try :
            data = [[trial, gid, node_type_id], cell_xyz, [el_dist, amp, spikes]]

            if include_sin:
                data = data + [[fq * 1.0]]

            table.loc[run_id] = list(itertools.chain.from_iterable(data))

        except:
            print run_id, [el_dist, amp, spikes]
            raise

    dt = conf['run']['dt']
    tstop = conf['run']['tstop']
    ex_waveform = stimx_waveform_factory(conf['inputs']['Extracellular_Stim']['waveform'])
    ex_delay = ex_waveform.delay
    ex_dur = ex_waveform.duration

    if include_spike_phase:

        print "Computing spike_threshold_time"
        data = cvh5['v/data'].value.T
        spike_tt = {}
        count = 0
        for vm in data:
            gid = gids[count]
            # print gid
            spike_tt[gid] = extract_spike_threshold_t(vm, dt=dt, ex_delay=ex_delay, ex_dur=ex_dur, tstop=tstop)
            count += 1

        print "Computing spike_phase"
        data = cvh5['vext/data'].value.T
        spike_phase = {}
        count = 0
        for vext in data:
            gid = gids[count]
            spike_ndx = get_spike_ndx(spike_tt[gid], dt=dt, ex_delay=ex_delay, ex_dur=ex_dur)
            phase_var = hilbert_transform(vext, dt=dt, ex_delay=ex_delay, ex_dur=ex_dur)
            spike_phase[gid] = [phase_var[i] for i in spike_ndx]
            count = count + 1

        additinal_table = pd.DataFrame(columns=['spike_phase', 'spike_threshold_t', 'gid'])

        for gid in gids:
            additinal_table = additinal_table.append({ "gid": gid,
                                              "spike_threshold_t": spike_tt[gid],
                                              "spike_phase": spike_phase[gid]
                                              }, ignore_index=True)

    Final_table = pd.merge(additinal_table, table, on=['gid'], how='outer')
    filename = nu.get_table_filename(network_id, amp, trial, fq)
    print 'Data collected. Writing to {}...'.format(filename)
    fpath = nu.get_table_dir(input_type, stim_type, model_type, saved_data, filename)
    print "writing to:", fpath
    nu.write_table_h5(fpath, Final_table, attrs={'has_sin': include_sin,
                                                 'has_spike_phase_analysis': include_spike_phase})
    print 'Done.'


def extract_spike_threshold_t(row, **kwargs):
    'We use Allensdk spike feature extractor to find the spike threshold time'
    dt = kwargs['dt']
    tstop = kwargs['tstop']
    voltage = np.asarray(row)
    time = np.arange(0, tstop, dt)
    time = time / 1000.
    spike_threshold_t = []

    if np.isnan(voltage).any():
        print "Nan values in the trace"
    else:
        sweep = EphysSweepFeatureExtractor(t=time, v=voltage, start=0, end=time[-1])
        sweep.process_spikes()
        all_spike_features = sweep.spike_feature_keys()

        if len(all_spike_features) > 0:
            features_dict = {}
            for j in range(0, len(all_spike_features)):
                features_dict[all_spike_features[j]] = sweep.spike_feature(all_spike_features[j])
            spike_threshold_t = [i * 1000 for i in features_dict['threshold_t'].tolist()]

    return spike_threshold_t


def hilbert_transform(row, **kwargs):
    'We compute the hilbert transform for the whole period when extra_stim is applied. After making the table\
    then we can cut the first and last 2s. But the table is build for the whole period of ex_stim applied'
    dt = kwargs['dt']
    ex_delay = kwargs['ex_delay']
    ex_dur = kwargs['ex_dur']
    t_start = int((ex_delay) / dt)
    t_end = int((ex_delay + ex_dur) / dt)
    var = row[t_start:t_end]
    n_points = int(1. / (dt * (0.001)))  # Number of data points in 1second

    if np.isnan(var).any():
        print "Nan values for phase computation"
        return []
    else:
        analytic_var = hilbert(var)
        amplitude_envelope_var = np.abs(analytic_var)
        instantaneous_phase_var = np.unwrap(np.angle(analytic_var))
        phase_var = [(x / (2.0 * np.pi) - int(x / (2.0 * np.pi))) * 2.0 * np.pi - np.pi for x in instantaneous_phase_var]
        instantaneous_frequency_var = (np.diff(instantaneous_phase_var) / (2.0 * np.pi) * n_points)

        return phase_var


def get_spike_ndx(spike_threshold_time, **kwargs):
        dt = kwargs['dt']
        ex_delay = kwargs['ex_delay']
        ex_dur = kwargs['ex_dur']
        t_start = ex_delay
        t_end = ex_delay + ex_dur
        time = np.arange(t_start, t_end, dt)

        ndx = []
        for t in spike_threshold_time:
            ndx = ndx + np.where((time > (t - 0.000001)) & (time < (t + 0.000001)))[0].tolist()

        return ndx






