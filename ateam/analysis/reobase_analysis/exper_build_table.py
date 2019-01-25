######################################################
# Authors: Fahimeh Baftizadeh
# Date created: 9/1/2017
######################################################

import glob
import numpy as np
import pandas as pd
import itertools
import allensdk
import reobase_analysis.reobase_utils as ru
import reobase_analysis.sin_utils as su
from scipy.signal import hilbert, chirp
from reobase_analysis.reobase_utils import StimType, ModelType, InputType
from reobase_analysis.stimxwaveform import stimx_waveform_factory
from reobase_analysis.stimxwaveform import iclamp_waveform_factory
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor

"""
Script to build h5 files containing organized run info for all runs of a particular amplitude
Will include all electrodes. Has no overwrite protection.

Will attempt to include vm information if the stim type is DC. This is only done for DC because it only makes sense for
a constant input without external input. This allows us to look at the subthreshold response.
If the cell is active (> 1 spikes) then the post-stim vm value is NaN
"""

def bandpass_filter():

    pass



def build_sin_dc(cell_gid, input_type, stim_type, model_type, inputs, trial, include_delta_vm=True,
                 include_sin=True, include_vm_phase_analysis=True, include_vext_phase_analysis=True,
                 include_vi_phase_analysis=True, include_spike_phase=True, saved_data=False,
                 control_simulation = True):
    cell_csv_pattern = '/*_cel[ls]*csv'
    cell_out_dir = ru.get_output_dir(input_type, stim_type, model_type, cell_gid, saved_data)
    include_iclamp = input_type == InputType.EXTRASTIM_INTRASTIM
    additional_cols = ru.resolve_additional_cols(include_delta_vm, include_sin, include_iclamp,
                                                 include_vm_phase_analysis, include_vext_phase_analysis,
                                                 include_vi_phase_analysis, include_spike_phase)
    print additional_cols
    print "Iclamp and stimx have the same delay and dur"
    for amp in inputs:
        print 'Build data table for ex_amp = {}...'.format(amp)
        all_dir_name = ru.resolve_output_aggregates(include_delta_vm, include_sin, include_iclamp, amp, trial)
        amp_output =  ru.concat_path(cell_out_dir, all_dir_name)
        table = ru.build_df(additional_cols)
        vm_rest = np.NaN
        for out_dir in glob.iglob(amp_output):
                dir_name = out_dir.split('/')[-1]
                el = int([x for x in dir_name.split('_') if x.startswith('el')][0][2:])
                fq = int([x for x in dir_name.split('_') if x.startswith('freq')][0][4:])
                ic_amp = int([x for x in dir_name.split('_') if x.startswith('icamp')][0][5:]) * 0.001 if "icamp" in dir_name else None

                config_path = ru.get_config_resolved_path(stim_type, out_dir, el, amp, fq, ic_amp)
                run_id = ru.resolve_run_id(cell_gid, el, amp, fq, ic_amp)
                cvh5_for_vext_trace = ru.get_cv_files(out_dir, [0])[0]

                conf = ru.get_json_from_file(config_path)
                electrodes_dir = ru.get_dir_root(saved_data) + '/' + '/'.join(conf['output']['electrodes_dir'].split('/')[2:])
                # electrodes_dir = conf['output']['electrodes_dir']

                in_waveform = iclamp_waveform_factory(conf)
                ex_waveform = stimx_waveform_factory(conf)
                ex_delay = ex_waveform.delay
                ex_dur = ex_waveform.duration
                if stim_type == StimType.SIN:
                    in_delay = None
                    in_dur = None
                elif stim_type == StimType.SIN_DC:
                    in_delay = in_waveform.delay
                    in_dur = in_waveform.duration

                if control_simulation:
                    control_outdir = ru.get_intrastim_output_dir(cell_gid, "intrastim", "dc", model_type, in_waveform.amp, trial, saved_data)
                    control_cvh5 = ru.get_cv_files(control_outdir, [0])[0]
                    cvh5_for_vm_trace = control_cvh5
                else:
                    cvh5_for_vm_trace = cvh5_for_vext_trace

                spikes = cvh5_for_vm_trace['spikes'].value

                el_xyz = ru.get_electrode_xyz(ru.get_electrode_path(electrodes_dir, cell_gid, el)).flatten()
                ru.get_electrode_path(electrodes_dir, cell_gid, el)
                cell_xyz = ru.get_cell_xyz(glob.glob(out_dir + cell_csv_pattern)[0]).flatten()
                el_dist = round(np.linalg.norm(el_xyz - cell_xyz))

                try:
                    data = [[trial, el], el_xyz, [el_dist, amp, spikes]]

                    if include_delta_vm:
                        vm_data = extract_vm_data(cvh5_for_vm_trace, stim_type, ex_delay, ex_dur, iclamp_delay = in_delay, iclamp_dur = in_dur) if include_delta_vm else []
                        vm_rest = vm_data.pop(0)
                        data = data + [vm_data]

                    if include_sin:
                        data = data + [[fq * 1.0]]

                    if include_iclamp:
                        data = data + [[ic_amp]]

                    if include_vm_phase_analysis:
                        vm_phase_analysis_data = extract_v_phase_analysis('vm', cvh5_for_vm_trace, cvh5_for_vext_trace, ex_delay, ex_dur, fq)
                        data = data + [vm_phase_analysis_data]

                    if include_vext_phase_analysis:
                        vext_phase_analysis_data = extract_v_phase_analysis('vext', cvh5_for_vm_trace, cvh5_for_vext_trace, ex_delay, ex_dur, fq)
                        data = data + [vext_phase_analysis_data]

                    if include_vi_phase_analysis:
                        vi_phase_analysis_data = extract_v_phase_analysis('vi', cvh5_for_vm_trace, cvh5_for_vext_trace, ex_delay, ex_dur, fq)
                        data = data + [vi_phase_analysis_data]

                    if include_spike_phase:
                        spike_threshold_data = extract_spike_threshold_t(cvh5_for_vm_trace)
                        spike_phase_data = extract_spike_phase('vext', cvh5_for_vm_trace, cvh5_for_vext_trace, ex_delay, ex_dur)
                        spike_threshold_phase_data = [spike_threshold_data] + [spike_phase_data]
                        data = data + [spike_threshold_phase_data]

                    table.loc[run_id] = list(itertools.chain.from_iterable(data))
                except:
                    print run_id, [el_dist, amp, spikes]
                    raise

        index_close_els = ru.get_index_close_els(cell_gid, input_type, stim_type, model_type, trial, saved_data)
        # table = table.drop(index_close_els)
        if include_vm_phase_analysis:
            table.loc[table['vm_phase'] < 0, 'vm_phase'] = table['vm_phase'] + 360
            table.loc[table["vi_phase"] < 0, "vi_phase"] = table["vi_phase"] + 360
            table.loc[table["vext_phase"] < 0, "vext_phase"] = table["vext_phase"] + 360
        if control_simulation:
            filename = ru.get_control_table_filename(cell_gid, amp,trial)
        else:
            filename = ru.get_table_filename(cell_gid, amp, trial)
        print 'Data collected. Writing to {}...'.format(filename)
        fpath = ru.get_table_dir(input_type,stim_type, model_type, filename)
        print "writing to:" , fpath
        ru.write_table_h5(fpath, table, attrs={'has_vm_data': include_delta_vm,
                                               'vm_rest': vm_rest,
                                               'has_iclamp': include_iclamp,
                                               'has_sin': include_sin,
                                               'has_vm_phase_analysis':include_vm_phase_analysis,
                                               'has_vext_phase_analysis': include_vext_phase_analysis,
                                               'has_vi_phase_analysis': include_vi_phase_analysis,
                                               'has_spike_phase_analysis': include_spike_phase})

        print 'Done.'


def extract_spike_threshold_t(cvh5_for_vm_trace):
    'We use Allensdk spike feature extractor to find the spike threshold time'
    spikes = cvh5_for_vm_trace['spikes'].value

    if len(spikes) < 1:
        spike_threshold_t = []
    else:
        voltage = cvh5_for_vm_trace['vm'].value
        dt = cvh5_for_vm_trace.attrs['dt']
        tstop = cvh5_for_vm_trace.attrs['tstop']
        time = np.arange(0,tstop,dt)
        time = time / 1000.
        sweep = EphysSweepFeatureExtractor(t=time, v=voltage, start=0, end=time[-1])
        sweep.process_spikes()
        all_spike_features = sweep.spike_feature_keys()

        features_dict = {}
        for j in range(0, len(all_spike_features)):
            features_dict[all_spike_features[j]] = sweep.spike_feature(all_spike_features[j])

        spike_threshold_t = features_dict['threshold_t']*1000.

    return spike_threshold_t


def hilbert_transform(varname, cvh5_for_vext_trace, ex_delay, ex_dur):
    'We compute the hilbert transform for the whole period when extra_stim is applied. After making the table\
    then we can cut the first and last 2s. But the table is build for the whole period of ex_stim applied'
    dt = cvh5_for_vext_trace.attrs['dt']
    var = cvh5_for_vext_trace[varname].value
    t_start = int((ex_delay) / dt)
    t_end = int((ex_delay + ex_dur) / dt)
    var = var[t_start:t_end]
    n_points = int(1. / (dt * (0.001))) #Number of data points in 1second

    analytic_var = hilbert(var)
    amplitude_envelope_var = np.abs(analytic_var)
    instantaneous_phase_var = np.unwrap(np.angle(analytic_var))
    phase_var = [(x / (2.0 * np.pi) - int(x / (2.0 * np.pi))) * 2.0 * np.pi - np.pi for x in instantaneous_phase_var]
    instantaneous_frequency_var = (np.diff(instantaneous_phase_var) / (2.0 * np.pi) * n_points)

    return phase_var, amplitude_envelope_var, instantaneous_frequency_var


def extract_spike_phase(varname, cvh5_for_vm_trace, cvh5_for_vext_trace, ex_delay, ex_dur):
    'We compute the spike phase for all the spikes from the beggining to end. After we build the table, then we can decide\
    which time window we keep for analyis. But as of now, we copute the spike phase for all the spikes'
    spikes = cvh5_for_vm_trace['spikes'].value

    if len(spikes) < 1:
        spike_phase = []
    else:
        dt = cvh5_for_vm_trace.attrs['dt']
        t_start = ex_delay
        t_end = ex_delay + ex_dur
        time = np.arange(t_start, t_end, dt)

        spike_threshold_time = extract_spike_threshold_t(cvh5_for_vm_trace)
        phase_var, amplitude_envelope_var, instantaneous_frequency_var = hilbert_transform(varname, cvh5_for_vext_trace, ex_delay, ex_dur)

        df = pd.DataFrame(columns=["spike_phase", 'time'])
        df['spike_phase'] = phase_var
        df["time"] = time

        ndx = []
        for t in spike_threshold_time:
            ndx = ndx + (df[(df['time'] > (t - 0.000001)) & (df['time'] < (t + 0.000001))].index.tolist())

        spike_phase = [df.loc[index]['spike_phase'] for index in ndx]

    return np.asarray(spike_phase)



def extract_v_phase_analysis(var_name ,cvh5_for_vm_trace, cvh5_for_vext_trace, ex_delay, ex_dur, freq):
    'This is subthreshold only analyis. When the cell is spiking this value is equal to NaN\
    We cut the first 4s and the final 2s for fitting the sinusoid'
    spikes = cvh5_for_vm_trace['spikes'].value

    if len(spikes) < 1 or var_name == 'vext':
        var_amp, var_phase, var_mean = su.fit_sin_h5(var_name, cvh5_for_vext_trace, ex_delay, ex_dur, freq)
    else:
        var_amp = np.NaN
        var_phase = np.NaN

    return [var_amp, var_phase]


def extract_vm_data(cvh5_for_vm_trace, stim_type, ex_delay, ex_dur, **kwargs):
    'This is a subthreshold feature only and when the cell is spiking, this value is equal to NaN\
    we cut the first 1s before the extracellular stim. and we copmute the v_rest as the final 500ms before extra_stim is applied\
    delta_vm is NaN is the cell is spiking, to compute the delta_vm for when the cell is not spiking if there is Iclamp\
    then we compute it 1s after Iclamp is injected. if there is no iclamp, then we compute delta_vm 1s after extra_stim is applied.'
    dt = cvh5_for_vm_trace.attrs['dt']
    vm = cvh5_for_vm_trace['vm'].value
    spikes = cvh5_for_vm_trace['spikes'].value
    get_step = lambda x: int(x / dt)

    vm_rest = vm[get_step(ex_delay - 500):get_step(ex_delay - 5)].mean()

    if len(spikes) < 1:  # ensure subthreshold or edge case
        if (stim_type == StimType.SIN):
            vm_stim = vm[get_step(ex_delay + 1000):get_step(ex_delay + ex_dur - 5)].mean()
            if (ex_dur) < 2000:
                print "ERROR1 in calculating vm_stim"

        elif stim_type == StimType.SIN_DC:
            in_delay = kwargs['iclamp_delay']
            in_dur = kwargs['iclamp_dur']
            vm_stim = vm[get_step(in_delay + 1000):get_step(in_delay + in_dur - 5)].mean()
            if (in_dur) < 2000:
                print "ERROR2 in calculating vm_stim"
        else:
            print "This STIMP_TYPE is not defined"
    else:
        vm_stim = np.NaN

    return [vm_rest, vm_stim, (vm_stim - vm_rest)]


def build_dc(cell_gid, ex_inputs, input_type, stim_type, model_type, trial):
    cell_csv_pattern = '/*_cel[ls]*csv'  # ridiculous pattern matching for old files called 1_cell.csv vs new ones called [gid].csv
    cell_out_dir = ru.get_output_dir(stim_type, model_type, cell_gid)
    include_delta_vm = stim_type == StimType.DC.value

    for amp in ex_inputs:
        print 'Build data table for amp = {}...'.format(amp)
        table = ru.build_df(ru.vm_cols if include_delta_vm else [])
        amp_output = ru.concat_path(cell_out_dir, ru.get_dc_dir_name(9999, amp, trial)).replace('9999', '*')
        vm_rest = np.NaN

        for out_dir in glob.iglob(amp_output):
            dir_name = out_dir.split('/')[-1]
            el = int([x for x in dir_name.split('_') if x.startswith('el')][0][2:]) # parse el from folder name..
            config_path = ru.get_config_resolved_path(stim_type, out_dir, el, amp)
            cvh5 = ru.get_cv_files(out_dir, [0])[0]
            conf = ru.get_json_from_file(config_path)
            electrodes_dir = ru.get_dir_root() + '/' + '/'.join( conf['output']['electrodes_dir'].split('/')[2:] )
            waveform = stimx_waveform_factory(conf)

            spikes = cvh5['spikes'].value
            el_xyz = ru.get_electrode_xyz(ru.get_electrode_path(electrodes_dir, cell_gid, el)).flatten()
            cell_xyz = ru.get_cell_xyz(glob.glob(out_dir + cell_csv_pattern)[0]).flatten()
            el_dist = round(np.linalg.norm(el_xyz - cell_xyz))
            run_id = ru.resolve_dc_run_id(cell_gid, el, amp)


            try:
                vm_data = extract_vm_data(cvh5, waveform.delay, waveform.duration) if include_delta_vm else []
                vm_rest = vm_data.pop(0)
                data = [[trial, el], el_xyz, [el_dist, amp, spikes], vm_data]
                table.loc[run_id] = list(itertools.chain.from_iterable(data))
            except:
                print run_id, [el_dist, amp, spikes]
                raise

        filename = ru.get_table_filename(cell_gid, amp, trial)
        print 'Data collected. Writing to {}...'.format(filename)
        fpath = ru.get_table_dir(stim_type, model_type, filename)
        ru.write_table_h5(fpath, table, attrs={'has_vm_data':include_delta_vm,'vm_rest': vm_rest})
        print 'Done.'
