######################################################
# Authors: Fahimeh Baftizadeh
# Date created: 9/1/2017
######################################################

import os, sys
import math
import glob
import json
import h5py as h5
import numpy as np
import pandas as pd
from enum import Enum
from scipy.optimize import leastsq
from scipy import fftpack
import ateam.analysis.brainstim.reobase_utils as ru
from .stimxwaveform import stimx_waveform_factory
from .stimxwaveform import iclamp_waveform_factory


def compute_fft(sig, time_step):
    N = len(sig)
    time_vec = np.arange(N)
    sig_fft = fftpack.fft(sig)
    power = np.abs(sig_fft)
    sample_freq = fftpack.fftfreq(sig.size, d=time_step)
    return sig_fft, sample_freq


def bandpass_filter(signal,lowcut_freq, highcut_freq, time_step):
    N = len(signal)
    time_vec = np.arange(N)
    sig_fft, sample_freq = compute_fft(signal, time_step)
    band_freq_fft = sig_fft.copy()
    band_freq_fft[(np.abs(sample_freq) > highcut_freq)] = 0
    band_freq_fft[(np.abs(sample_freq) < lowcut_freq)] = 0
    filtered_sig = fftpack.ifft(band_freq_fft)
    return np.real(filtered_sig)

def get_var(varname, cvh5):
    if varname == "vi":
        var = np.add(cvh5['vm'].value, cvh5['vext'].value)
    else:
        var = cvh5[varname].value
    return var

def get_dt(cvh5):
    return cvh5.attrs['dt']

def get_mean_trace(var_trace, ex_delay, ex_dur, freq, dt):
    # print "here3", freq
    var = var_trace
    t_start_analysis = ex_delay + 3000.
    t_end_analysis = ex_delay + ex_dur - 2000

    n_period = 1000. / freq  #Time for each cycle
    n_traces = int((t_end_analysis - t_start_analysis) / round(n_period))
    traces = []
    slice_tsteps = [int(i / dt) for i in np.arange(t_start_analysis, t_end_analysis+ 1, n_period)]
    for beg_trace, end_trace in zip(slice_tsteps[:-1], slice_tsteps[1:]):
        mean = np.mean(var[beg_trace:end_trace])
        traces.append(var[beg_trace:end_trace] - mean)

    mean_list_of_list = lambda a: sum(a) * 1. / len(a)
    mean_trace = map(mean_list_of_list, zip(*traces))
    return  mean_trace

def fit_sin(var_trace, ex_delay, ex_dur, freq, dt):
    # print "here2"
    data = get_mean_trace(var_trace, ex_delay, ex_dur, freq, dt)
    N = len(data)  # Number of data points
    t = np.linspace(0, 2 * np.pi, N)

    guess_mean = 0
    guess_std = 0.002
    guess_phase = 1

    # data_first_guess = guess_std * np.sin(t + guess_phase) + guess_mean
    optimize_func = lambda x: x[0] * np.sin(t + x[1]) + x[2] - data
    est_std, est_phase, est_mean = leastsq(optimize_func, [guess_std, guess_phase, guess_mean])[0]
    data_fit = est_std * np.sin(t + est_phase) + est_mean

    if est_std < 0:
        est_phase = est_phase + np.pi
        est_std *= -1.0

    est_phase = est_phase / (2 * np.pi) - int(est_phase / (2 * np.pi))
    est_phase = est_phase * 360

    if abs(est_phase) > 360:
        print "Estimated Phase can not be more than 360"
    # elif est_phase < -50:
    #     est_phase = est_phase + 360

    return est_std, est_phase, est_mean
    # return t, data

def get_fitted_sin(var_trace, ex_delay, ex_dur, freq, dt):
    # print "here2"
    data = get_mean_trace(var_trace, ex_delay, ex_dur, freq, dt)
    N = len(data)  # Number of data points
    t = np.linspace(0, 2 * np.pi, N)

    guess_mean = 0
    guess_std = 0.002
    guess_phase = 1

    # data_first_guess = guess_std * np.sin(t + guess_phase) + guess_mean
    optimize_func = lambda x: x[0] * np.sin(t + x[1]) + x[2] - data
    est_std, est_phase, est_mean = leastsq(optimize_func, [guess_std, guess_phase, guess_mean])[0]
    data_fit = est_std * np.sin(t + est_phase) + est_mean

    if est_std < 0:
        est_phase = est_phase + np.pi
        est_std *= -1.0

    est_phase = est_phase / (2 * np.pi) - int(est_phase / (2 * np.pi))
    est_phase = est_phase * 360

    if abs(est_phase) > 360:
        print "Estimated Phase can not be more than 360"
    # elif est_phase < -50:
    #     est_phase = est_phase + 360

    return data_fit



def get_mean_trace_h5(var_name ,cvh5, ex_delay, ex_dur, freq):
    'This is used for subthreshold analysis only. We cut the first 4s (which is ex_delay+3000) also the final 2s\
    which is ex_delay + ex_dur - 2000'
    dt = cvh5.attrs['dt']
    if var_name == 'vi':
        var = np.add(cvh5['vm'].value, cvh5['vext'].value)
    else:
        var = cvh5[var_name].value
    t_start_analysis = ex_delay + 3000.
    t_end_analysis = ex_delay + ex_dur - 2000

    n_period = 1000. / freq  #Time for each cycle
    n_traces = int((t_end_analysis - t_start_analysis) / round(n_period))
    traces = []
    slice_tsteps = [int(i / dt) for i in np.arange(t_start_analysis, t_end_analysis+ 1, n_period)]
    for beg_trace, end_trace in zip(slice_tsteps[:-1], slice_tsteps[1:]):
        mean = np.mean(var[beg_trace:end_trace])
        traces.append(var[beg_trace:end_trace] - mean)

    mean_list_of_list = lambda a: sum(a) * 1. / len(a)
    mean_trace = map(mean_list_of_list, zip(*traces))

    return  mean_trace


def fit_sin_h5(var_name ,cvh5, ex_delay, ex_dur, freq):
    data = get_mean_trace_h5(var_name ,cvh5, ex_delay, ex_dur, freq)
    N = len(data)  # Number of data points
    t = np.linspace(0, 2 * np.pi, N)

    guess_mean = 0
    guess_std = 0.002
    guess_phase = 1

    # data_first_guess = guess_std * np.sin(t + guess_phase) + guess_mean
    optimize_func = lambda x: x[0] * np.sin(t + x[1]) + x[2] - data
    est_std, est_phase, est_mean = leastsq(optimize_func, [guess_std, guess_phase, guess_mean])[0]
    data_fit = est_std * np.sin(t + est_phase) + est_mean

    if est_std < 0:
        est_phase = est_phase + np.pi
        est_std *= -1.0

    est_phase = est_phase / (2 * np.pi) - int(est_phase / (2 * np.pi))
    est_phase = est_phase * 360

    if abs(est_phase) > 360:
        print "Estimated Phase can not be more than 360"
    # elif est_phase < -50:
    #     est_phase = est_phase + 360

    return est_std, est_phase, est_mean


# def get_mean_trace(var_name, cell_id,input_type, stim_type, model_type, trial, el, amp, freq, ic_amp=None):
#
#     output_dir = ru.get_output_dir(input_type, stim_type, model_type, cell_id)
#     output_file = ru.get_dir_name(el, amp, freq, ic_amp, trial)
#     out_dir = ru.concat_path(output_dir, output_file)
#     cvh5 = ru.get_cv_files(out_dir, [0])[0]
#     config_path = ru.get_config_resolved_path(stim_type, out_dir, el, amp, freq, ic_amp)
#     conf = ru.get_json_from_file(config_path)
#     in_waveform = iclamp_waveform_factory(conf)
#
#     dt = cvh5.attrs['dt']
#     var = cvh5[var_name].value
#     t_start_analysis = in_waveform.delay + 1000.
#     t_end_analysis = in_waveform.delay + in_waveform.duration
#
#     n_period = 1000. / freq  #Time for each cycle
#     n_traces = int((t_end_analysis - t_start_analysis) / round(n_period))
#     traces = []
#     slice_tsteps = [int(i / dt) for i in np.arange(t_start_analysis, t_end_analysis+ 1, n_period)]
#     for beg_trace, end_trace in zip(slice_tsteps[:-1], slice_tsteps[1:]):
#         mean = np.mean(var[beg_trace:end_trace])
#         traces.append(var[beg_trace:end_trace] - mean)
#
#     mean_list_of_list = lambda a: sum(a) * 1. / len(a)
#     mean_trace = map(mean_list_of_list, zip(*traces))
#
#     return  mean_trace
#



# def fit_sin(var_name, cell_id,input_type, stim_type, model_type, trial, el, amp, freq, ic_amp=None):
#
#     data = get_mean_trace(var_name, cell_id,input_type, stim_type, model_type, trial, el, amp, freq, ic_amp)
#     N = len(data) #Number of data points
#     t = np.linspace(0, 2 * np.pi, N)
#
#     guess_mean = 0
#     guess_std = 0.002
#     guess_phase = 1
#
#     data_first_guess = guess_std * np.sin(t + guess_phase) + guess_mean
#     optimize_func = lambda x: x[0] * np.sin(t + x[1]) + x[2] - data
#     est_std, est_phase, est_mean = leastsq(optimize_func, [guess_std, guess_phase, guess_mean])[0]
#     data_fit = est_std * np.sin(t + est_phase) + est_mean
#
#     if est_std < 0:
#         est_phase = est_phase + np.pi
#         est_std *= -1.0
#
#     est_phase = est_phase / (2 * np.pi) - int(est_phase / (2 * np.pi))
#     est_phase = est_phase * 360
#
#     if abs(est_phase) > 360:
#         print "Estimated Phase can not be more than 360"
#     # elif est_phase < -50:
#     #     est_phase = est_phase + 360
#
#     return est_std, est_phase, est_mean