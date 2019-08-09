######################################################
# Authors: Fahimeh Baftizadeh, Taylor Connington
# Date created: 9/1/2017
######################################################

import glob
import numpy as np
import pandas as pd
import itertools
import allensdk
import ateam.analysis.brainstim.reobase_utils as ru
import ateam.analysis.brainstim.exper_utils as xu
import ateam.analysis.brainstim.config as config
import ateam.analysis.brainstim.sin_utils as su
from scipy.signal import hilbert, chirp
from .exper_utils import ElectrodeType
from .reobase_utils import StimType, ModelType, InputType
from .stimxwaveform import stimx_waveform_factory
from .stimxwaveform import iclamp_waveform_factory
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor
from neuroanalysis.miesnwb import MiesNwb

"""
Script to build h5 files containing organized run info for all runs of a particular amplitude
Will include all electrodes. Has no overwrite protection.
"""

##########Building experimental table###############

# Each table contains stimulus data, vext_phase
# analysis, and all spike analysis, if you dont
# want any of these parts simply set it to false
# on the command line

def build_expr_table(exp_id, sampling_freq, lowcut_freq, highcut_freq, saved_data = False, include_ex_stimulus_data = True,
                     include_in_stimulus_data = True, include_vext_phase_analysis = True,
                     include_vi_phase_analysis = True, include_spike_phase = True, include_vm_phase_analysis= True):

    experimental_cols = xu.resolve_additional_cols(include_ex_stimulus_data, include_in_stimulus_data,
                                                   include_vext_phase_analysis, include_vi_phase_analysis,
                                                   include_spike_phase, include_vm_phase_analysis)
    print "Building table for experiment:{}_{}".format(exp_id, sampling_freq)
    ####### First we generate an empty data frame with all the expected col names #######
    table = xu.build_df(experimental_cols)
    ####### We read all the necessary info from the config file #######
    original_config = xu.get_config_path(exp_id, sampling_freq, saved_data)
    config.print_resolved(config.build(original_config))
    resolved_config_path = xu.get_config_resolved_path(exp_id, sampling_freq, saved_data)
    conf = config.build(resolved_config_path)

    ex_el = xu.get_el_prop(conf, ElectrodeType.EXTRA)
    if (len(ex_el['id']) != len(ex_el['distance'])):
        raise ValueError('For each ex electrode, you should provide a distance')

    in_el = xu.get_el_prop(conf, ElectrodeType.INTRA)
    in_el_id = in_el['id']

    stim_el = xu.get_el_prop(conf, ElectrodeType.STIM)
    stim_el_id = stim_el['id']

    nwb_file_name = xu.get_nwb_filename(exp_id, sampling_freq)
    nwb_file = "/".join([xu.get_nwb_dir(conf), nwb_file_name])
    nwb = MiesNwb(nwb_file)

    sweeps = conf['sweep_numbers']
    control_sweeps = conf['control_sweep_numbers']
    number_of_sweeps_after_control = conf['number_of_sweeps_after_control']

    for sweep in sweeps:
        subthreshold = True
        control = False
        ####### extract the extracellular and intracellular injection properties #######
        in_amp, in_frequency, in_delay, in_dur = extract_stimulus_prop(nwb, conf, sweep, in_el_id, stimulus_type="in")
        ex_amp, ex_frequency, ex_delay, ex_dur =  extract_stimulus_prop(nwb, conf, sweep, stim_el_id, stimulus_type = "ex")

        if ex_amp == 0 :
            control = True
            if not sweep in control_sweeps: raise ValueError('it is a control sweep but its number is not inside the control list. Check sweep: {}'.format(sweep))
            related_sweeps = get_related_sweeps_to_control(sweep, sweeps, number_of_sweeps_after_control)
            print "For control sweep:", sweep, "these are related sweeps with different frequencies:", related_sweeps
        else:
            related_sweeps = [sweep]

        total_sweep_time = nwb.contents[sweep][stim_el_id]['primary'].duration
        total_sweep_length = total_sweep_time /(1./(sampling_freq * 1000))
        sweep_time_array = extract_timeseries_of_experiment(total_sweep_time, sampling_freq) # in sec
        dt = (sweep_time_array[1] - sweep_time_array[0]) * 1000 #in ms
        ex_delay = ex_delay * 1000 #in ms
        ex_dur = ex_dur * 1000 # in ms
        in_delay = in_delay * 1000 #in ms
        in_dur = in_dur * 1000 # in ms
        vi_trace = nwb.contents[sweep][in_el_id]['primary'].data * 1000. # in mV
        if (len(vi_trace) != total_sweep_length):
            raise ValueError('vi_trace should have a lenght of {}'.format(total_sweep_length))
        spike_tt = extract_spike_threshold_t_expr(sweep_time_array, vi_trace)
        spike_tt= [stt for stt in spike_tt if (stt >= ex_delay) & (stt <= ex_delay + ex_dur)] # Discard all the spikes before and after extracellular injection


        # # Decide to perform sub or supra analysis
        if len(spike_tt) > 0:
            subthreshold = False

        i = 0
        for exel in ex_el['id']:
            distance = ex_el['distance'][i]
            for s in related_sweeps:
                sweep_for_vext = s
                signal = nwb.contents[sweep_for_vext][exel]['primary'].data * 1000 # in mV
                vext_trace = su.bandpass_filter(signal, lowcut_freq, highcut_freq, dt / 1000.) # dt should be in sec
                run_id = xu.resolve_run_id(sweep, exel, ex_amp, ex_frequency, in_amp)

                if sweep != s: # We just update the run_id and frequency
                    new_ex_amp, ex_frequency, ex_delay, ex_dur = extract_stimulus_prop(nwb, conf, s, stim_el_id,stimulus_type="ex")
                    ex_delay = ex_delay * 1000  # in ms
                    ex_dur = ex_dur * 1000  # in ms
                    run_id = xu.resolve_run_id(sweep, exel, ex_amp, ex_frequency, in_amp)
                    run_id = str(s) + "_" + run_id

                if (len(vext_trace) != total_sweep_length):
                    raise ValueError('vext_trace should have a lenght of {}'.format(total_sweep_length))


                try:
                    data = write_data_rows(sweep, exel, stim_el_id, in_el_id, distance, spike_tt, ex_amp, ex_frequency,
                                       ex_dur, ex_delay , in_amp,in_dur, in_delay,subthreshold, control,vext_trace,dt,
                                       vi_trace, include_ex_stimulus_data, include_in_stimulus_data, include_vext_phase_analysis,
                                       include_vi_phase_analysis, include_spike_phase, include_vm_phase_analysis)


                    table.loc[run_id] = list(itertools.chain.from_iterable(data))

                except:
                    print run_id
                    raise
            i += 1

    if include_vi_phase_analysis:
            table.loc[table["vi_phase"] < 0, "vi_phase"] = table["vi_phase"] + 360
    if include_vext_phase_analysis:
            table.loc[table["vext_phase"] < 0, "vext_phase"] = table["vext_phase"] + 360
    if include_vm_phase_analysis:
            table.loc[table["vm_phase"] < 0, "vm_phase"] = table["vm_phase"] + 360

    filename = xu.get_table_filename(exp_id, sampling_freq)
    print 'Data collected. Writing to {}...'.format(filename)
    fpath = xu.get_table_dir(filename, saved_data)
    print "writing to:" , fpath
    xu.write_table_h5(fpath, table, attrs={'has_ex_stimulus_data':include_ex_stimulus_data,
                                            'has_in_stimulus_data': include_in_stimulus_data,
                                            'include_vext_phase_analysis': include_vext_phase_analysis,
                                            'include_vi_phase_analysis':include_vi_phase_analysis,
                                            'include_vm_phase_analysis':include_vm_phase_analysis,
                                            'include_spike_phase':include_spike_phase})

    print 'Done.'


def write_data_rows(sweep, exel, stim_el_id, in_el_id, distance, spike_tt, ex_amp, ex_frequency, ex_dur, ex_delay , in_amp,
                    in_dur, in_delay,subthreshold, control,vext_trace,dt, vi_trace, include_ex_stimulus_data, include_in_stimulus_data,
                    include_vext_phase_analysis, include_vi_phase_analysis, include_spike_phase, include_vm_phase_analysis):

    data = [[sweep, exel, stim_el_id, in_el_id, distance], [spike_tt]]
    if include_ex_stimulus_data:
        data = data + [[ex_amp, ex_frequency, ex_dur, ex_delay]]
    if include_in_stimulus_data:
        data = data + [[in_amp, in_dur, in_delay]]

    if include_vext_phase_analysis:
        vext_phase_analysis_data = extract_v_phase_analysis_expr(subthreshold, control, 'vext', vext_trace, ex_delay,
                                                                 ex_dur, ex_frequency, dt)
        data = data + [vext_phase_analysis_data]

    if include_vi_phase_analysis:
        vi_phase_analysis_data = extract_v_phase_analysis_expr(subthreshold, control, 'vi', vi_trace, ex_delay, ex_dur,
                                                               ex_frequency, dt)
        data = data + [vi_phase_analysis_data]

    if include_spike_phase:
        spike_phase_data = extract_spike_phase_expr(subthreshold, spike_tt, vext_trace, ex_delay, ex_dur, dt)
        data = data + [spike_phase_data]

    if include_vm_phase_analysis:
        if subthreshold and distance == 50:
            vm_trace = vi_trace - vext_trace
            avg_vm = get_avg_v(vm_trace, in_delay, in_dur, dt)
            vm_phase_analysis_data = extract_v_phase_analysis_expr(subthreshold, control, 'vm', vm_trace, ex_delay, ex_dur,
                                                               ex_frequency, dt)
            vm_phase_analysis_data = vm_phase_analysis_data + [avg_vm]
        else:
            vm_phase_analysis_data = [np.NAN, np.NAN, np.NaN]
        data = data + [vm_phase_analysis_data]

    return  data

def get_avg_v(vm_trace, delay, dur, dt):
    t_start = int((delay) / dt)
    t_end = int((delay + dur) / dt)
    return np.mean(vm_trace[t_start:t_end])

def get_related_sweeps_to_control(control_sweep, sweeps, number_of_sweeps_after_control):
    ndx = sweeps.index(control_sweep)
    return [sweeps[i] for i in [ndx+j for j in range(1, number_of_sweeps_after_control+1)]]


def extract_timeseries_from_h5(cvh5, var):
    voltage = cvh5[var].value
    dt = cvh5.attrs['dt']
    tstop = cvh5.attrs['tstop']
    time = np.arange(0,tstop,dt)
    time = time / 1000.
    return time , voltage

def extract_timeseries_of_experiment(total_sweep_time, sampling_freq):
    dt = 1. / (sampling_freq * 1000)
    time = np.arange(0, total_sweep_time, dt)
    return time


def extract_v_phase_analysis_expr(subthreshold, control, var_name ,var_trace, ex_delay, ex_dur, freq, dt):
    'This is subthreshold only analyis. When the cell is spiking this value is equal to NaN\
    We cut the first 4s and the final 2s for fitting the sinusoid'
    # spikes = cvh5_for_vm_trace['spikes'].value
    if control:
        return [np.NaN, np.NaN]

    if subthreshold or var_name == 'vext':
        var_amp, var_phase, var_mean = su.fit_sin(var_trace, ex_delay, ex_dur, freq, dt)
    else:
        var_amp = np.NaN
        var_phase = np.NaN

    return [var_amp, var_phase]


def extract_spike_threshold_t_expr(time, var_timeseries):
    spike_threshold_t = []

    if np.isnan(var_timeseries).any():
        print "Nan values in the trace"
    else:
        sweep = EphysSweepFeatureExtractor(t=time, v=var_timeseries, start=0, end=time[-1])
        sweep.process_spikes()
        all_spike_features = sweep.spike_feature_keys()

        if len(all_spike_features) > 0:
            features_dict = {}
            for j in range(0, len(all_spike_features)):
                features_dict[all_spike_features[j]] = sweep.spike_feature(all_spike_features[j])
            spike_threshold_t = [i * 1000 for i in features_dict['threshold_t'].tolist()]

    return spike_threshold_t

def extract_stimulus_prop(nwb, conf, sweep, el_id, stimulus_type):
    stimulus_description = nwb.contents[sweep][el_id].stimulus.description
    BB = conf['stimulus description'][str(stimulus_description)]
    amp = nwb.contents[sweep][el_id].stimulus.items[BB].amplitude  # in nanoamp
    frequency = None

    if stimulus_type == "ex":
        amp = amp * 1e9
        if hasattr(nwb.contents[sweep][el_id].stimulus.items[BB], 'frequency'):
            frequency = nwb.contents[sweep][el_id].stimulus.items[BB].frequency
        else:
            frequency = 0
    if stimulus_type == "in":
        amp = amp * 1e12
        if hasattr(nwb.contents[sweep][el_id].stimulus.items[BB], 'frequency'):
            raise ValueError('Intracellular injection can not jave frequency attr')

    delay = nwb.contents[sweep][el_id].stimulus.items[BB].global_start_time
    dur = nwb.contents[sweep][el_id].stimulus.items[BB].duration

    return amp, frequency, delay, dur

def hilbert_transform_expr(var_trace, ex_delay, ex_dur,dt):
    'We compute the hilbert transform for the whole period when extra_stim is applied. After making the table\
    then we can cut the first and last 2s. But the table is build for the whole period of ex_stim applied'
    # dt = cvh5_for_vext_trace.attrs['dt']
    # var = cvh5_for_vext_trace[varname].value
    var = var_trace
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

def extract_spike_phase_expr(subthreshold, spike_tt, var_trace, ex_delay, ex_dur, dt):
    'We compute the spike phase for all the spikes from the beggining to end. After we build the table, then we can decide\
    which time window we keep for analyis. But as of now, we copute the spike phase for all the spikes'
    # spikes = cvh5_for_vm_trace['spikes'].value

    if subthreshold:
        spike_phase = []
    else:
        t_start = ex_delay
        t_end = ex_delay + ex_dur
        time = np.arange(t_start, t_end, dt)

        spike_threshold_time = spike_tt
        phase_var, amplitude_envelope_var, instantaneous_frequency_var = hilbert_transform_expr( var_trace, ex_delay, ex_dur, dt)

        df = pd.DataFrame(columns=["spike_phase", 'time'])
        df['spike_phase'] = phase_var
        df["time"] = time

        ndx = []
        for t in spike_threshold_time:
            ndx = ndx + (df[(df['time'] > (t - 0.000001)) & (df['time'] < (t + 0.000001))].index.tolist())

        spike_phase = [df.loc[index]['spike_phase'] for index in ndx]

    return [spike_phase]


##########Building simulation table###############

# Each table contains stimulus data, vext_phase
# analysis, and all spike analysis, if you dont
# want any of these parts simply set it to false
# on the command line


# def build_sin_dc(cell_gid, input_type, stim_type, model_type, inputs, trial, include_delta_vm=True,
#                  include_sin=True, include_vm_phase_analysis=True, include_vext_phase_analysis=True,
#                  include_vi_phase_analysis=True, include_spike_phase=True, saved_data=False,
#                  control_simulation = True):
#     cell_csv_pattern = '/*_cel[ls]*csv'
#     cell_out_dir = ru.get_output_dir(input_type, stim_type, model_type, cell_gid, saved_data)
#     include_iclamp = input_type == InputType.EXTRASTIM_INTRASTIM
#     additional_cols = ru.resolve_additional_cols(include_delta_vm, include_sin, include_iclamp,
#                                                  include_vm_phase_analysis, include_vext_phase_analysis,
#                                                  include_vi_phase_analysis, include_spike_phase)
#     print additional_cols
#     print "Iclamp and stimx have the same delay and dur"
#     for amp in inputs:
#         print 'Build data table for ex_amp = {}...'.format(amp)
#         all_dir_name = ru.resolve_output_aggregates(include_delta_vm, include_sin, include_iclamp, amp, trial)
#         amp_output =  ru.concat_path(cell_out_dir, all_dir_name)
#         table = ru.build_df(additional_cols)
#         vm_rest = np.NaN
#         for out_dir in glob.iglob(amp_output):
#                 dir_name = out_dir.split('/')[-1]
#                 el = int([x for x in dir_name.split('_') if x.startswith('el')][0][2:])
#                 fq = int([x for x in dir_name.split('_') if x.startswith('freq')][0][4:])
#                 ic_amp = int([x for x in dir_name.split('_') if x.startswith('icamp')][0][5:]) * 0.001 if "icamp" in dir_name else None

#                 config_path = ru.get_config_resolved_path(stim_type, out_dir, el, amp, fq, ic_amp)
#                 run_id = ru.resolve_run_id(cell_gid, el, amp, fq, ic_amp)
#                 cvh5_for_vext_trace = ru.get_cv_files(out_dir, [0])[0]

#                 conf = ru.get_json_from_file(config_path)
#                 electrodes_dir = ru.get_dir_root(saved_data) + '/' + '/'.join(conf['output']['electrodes_dir'].split('/')[2:])
#                 # electrodes_dir = conf['output']['electrodes_dir']

#                 in_waveform = iclamp_waveform_factory(conf)
#                 ex_waveform = stimx_waveform_factory(conf)
#                 ex_delay = ex_waveform.delay
#                 ex_dur = ex_waveform.duration
#                 if stim_type == StimType.SIN:
#                     in_delay = None
#                     in_dur = None
#                 elif stim_type == StimType.SIN_DC:
#                     in_delay = in_waveform.delay
#                     in_dur = in_waveform.duration

#                 if control_simulation:
#                     control_outdir = ru.get_intrastim_output_dir(cell_gid, "intrastim", "dc", model_type, in_waveform.amp, trial, saved_data)
#                     control_cvh5 = ru.get_cv_files(control_outdir, [0])[0]
#                     cvh5_for_vm_trace = control_cvh5
#                 else:
#                     cvh5_for_vm_trace = cvh5_for_vext_trace

#                 spikes = cvh5_for_vm_trace['spikes'].value

#                 el_xyz = ru.get_electrode_xyz(ru.get_electrode_path(electrodes_dir, cell_gid, el)).flatten()
#                 ru.get_electrode_path(electrodes_dir, cell_gid, el)
#                 cell_xyz = ru.get_cell_xyz(glob.glob(out_dir + cell_csv_pattern)[0]).flatten()
#                 el_dist = round(np.linalg.norm(el_xyz - cell_xyz))

#                 try:
#                     data = [[trial, el], el_xyz, [el_dist, amp, spikes]]

#                     if include_delta_vm:
#                         vm_data = extract_vm_data(cvh5_for_vm_trace, stim_type, ex_delay, ex_dur, iclamp_delay = in_delay, iclamp_dur = in_dur) if include_delta_vm else []
#                         vm_rest = vm_data.pop(0)
#                         data = data + [vm_data]

#                     if include_sin:
#                         data = data + [[fq * 1.0]]

#                     if include_iclamp:
#                         data = data + [[ic_amp]]

#                     if include_vm_phase_analysis:
#                         vm_phase_analysis_data = extract_v_phase_analysis('vm', cvh5_for_vm_trace, cvh5_for_vext_trace, ex_delay, ex_dur, fq)
#                         data = data + [vm_phase_analysis_data]

#                     if include_vext_phase_analysis:
#                         vext_phase_analysis_data = extract_v_phase_analysis('vext', cvh5_for_vm_trace, cvh5_for_vext_trace, ex_delay, ex_dur, fq)
#                         data = data + [vext_phase_analysis_data]

#                     if include_vi_phase_analysis:
#                         vi_phase_analysis_data = extract_v_phase_analysis('vi', cvh5_for_vm_trace, cvh5_for_vext_trace, ex_delay, ex_dur, fq)
#                         data = data + [vi_phase_analysis_data]

#                     if include_spike_phase:
#                         spike_threshold_data = extract_spike_threshold_t(cvh5_for_vm_trace)
#                         spike_phase_data = extract_spike_phase('vext', cvh5_for_vm_trace, cvh5_for_vext_trace, ex_delay, ex_dur)
#                         spike_threshold_phase_data = [spike_threshold_data] + [spike_phase_data]
#                         data = data + [spike_threshold_phase_data]

#                     table.loc[run_id] = list(itertools.chain.from_iterable(data))
#                 except:
#                     print run_id, [el_dist, amp, spikes]
#                     raise

#         index_close_els = ru.get_index_close_els(cell_gid, input_type, stim_type, model_type, trial, saved_data)
#         # table = table.drop(index_close_els)
#         if include_vm_phase_analysis:
#             table.loc[table['vm_phase'] < 0, 'vm_phase'] = table['vm_phase'] + 360
#             table.loc[table["vi_phase"] < 0, "vi_phase"] = table["vi_phase"] + 360
#             table.loc[table["vext_phase"] < 0, "vext_phase"] = table["vext_phase"] + 360
#         if control_simulation:
#             filename = ru.get_control_table_filename(cell_gid, amp,trial)
#         else:
#             filename = ru.get_table_filename(cell_gid, amp, trial)
#         print 'Data collected. Writing to {}...'.format(filename)
#         fpath = ru.get_table_dir(input_type,stim_type, model_type, filename)
#         print "writing to:" , fpath
#         ru.write_table_h5(fpath, table, attrs={'has_vm_data': include_delta_vm,
#                                                'vm_rest': vm_rest,
#                                                'has_iclamp': include_iclamp,
#                                                'has_sin': include_sin,
#                                                'has_vm_phase_analysis':include_vm_phase_analysis,
#                                                'has_vext_phase_analysis': include_vext_phase_analysis,
#                                                'has_vi_phase_analysis': include_vi_phase_analysis,
#                                                'has_spike_phase_analysis': include_spike_phase})

#         print 'Done.'




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

# def hilbert_transform(time, var_timeseries, ex_delay, ex_dur):
#     'We compute the hilbert transform for the whole period when extra_stim is applied. After making the table\
#     then we can cut the first and last 2s. But the table is build for the whole period of ex_stim applied'
#     dt = (time[1]- time[0]) / 1000.
#     t_start = int((ex_delay) / dt)
#     t_end = int((ex_delay + ex_dur) / dt)
#     var = var_timeseries[t_start:t_end]
#     n_points = int(1. / (dt * (0.001))) #Number of data points in 1second
#
#     analytic_var = hilbert(var)
#     amplitude_envelope_var = np.abs(analytic_var)
#     instantaneous_phase_var = np.unwrap(np.angle(analytic_var))
#     phase_var = [(x / (2.0 * np.pi) - int(x / (2.0 * np.pi))) * 2.0 * np.pi - np.pi for x in instantaneous_phase_var]
#     instantaneous_frequency_var = (np.diff(instantaneous_phase_var) / (2.0 * np.pi) * n_points)
#
#     return phase_var, amplitude_envelope_var, instantaneous_frequency_var


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


def compute_STA(cvh5, fq):
    dt = cvh5.attrs['dt']
    spike_threshold_time = extract_spike_threshold_t(cvh5)
    vext = cvh5['vext'].value

    # Windows from Costas' paper
    if fq == 1:
        window = 1
    elif fq == 8:
        window = 0.5
    elif fq == 30:
        window = 0.1
    else:
        window = 1./fq

    window = window * 1000
    traces = []
    for time in spike_threshold_time:
        beg_trace = (time - window) / dt
        end_trace = (time + window) / dt
        traces.append(vext[int(round(beg_trace)):int(round(end_trace))])
    traces = np.array(traces)
    mean_trace = np.mean(traces, axis=0)

    return mean_trace


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

