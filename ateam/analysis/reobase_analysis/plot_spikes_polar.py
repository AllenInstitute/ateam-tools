######################################################
# Authors: Fahimeh Baftizadeh, Taylor Connington
# Date created: 9/1/2017
######################################################

from scipy.signal import hilbert, chirp
import numpy as np
import matplotlib.pyplot as plt
import reobase_analysis.reobase_utils as ru
import reobase_analysis.tchelpers as tc
import reobase_analysis.plot_polar_angle as pol
import pycircstat
import pandas as pd
import random


def filter_list(row, **kwargs):
    t_start_analysis = 4000.0
    t_end_analysis = 12000.0
    varname1 = kwargs['varname1']
    varname2 = kwargs['varname2']
    ndx = np.where((row[varname1] >= t_start_analysis) & (row[varname1] <= t_end_analysis))
    l = [row[varname2][i] for i in ndx]
    flat_list = [item for sublist in l for item in sublist]
    return flat_list


def phase_correction(row):
    '''Shifts Hilbert transform result from -Pi to Pi to 0 to 2*Pi.'''
    phases_shifted = [x + (1.5 * np.pi) for x in row['spike_phase_window']]
    return [(x / (2.0 * np.pi) - int( x / (2.0 * np.pi))) * 2.0 * np.pi for x in phases_shifted]


def add_corrected_spikes_columns(table):
    '''Add columns to dataframe for analysis. Cuts out first 4 seconds and last 3 seconds for analysis.'''
    table['spikes_window'] = table.apply(filter_list, varname1="spikes", varname2="spikes", axis=1)
    table['spike_threshold_t_window'] = table.apply(filter_list, varname1="spikes", varname2="spike_threshold_t", axis=1)
    table['spike_phase_window'] = table.apply(filter_list, varname1="spikes", varname2="spike_phase", axis=1)
    table['spike_phase_correct_window'] = table.apply(phase_correction, axis=1)
    table['num_spikes_window'] = table.apply(lambda row: len(row['spikes_window']), axis=1)
    return table


def generate_spike_analysis_table(gid_list, input_type, stim_type, model_type, inputs, trial):
    '''Reads the tables for each cell ID and concatenates them. Adds additional columns
    for spike analysis.'''
    table = pd.DataFrame()   
    for cell_gid in gid_list:
        t = pol.generate_theta(amp=inputs, input_type=input_type, stim_type=stim_type, 
                               model_type=model_type, trial=trial, cell_gid=cell_gid)
        table = pd.concat([table, t])
    table = add_corrected_spikes_columns(table)
    return table


def plot_polar_histogram_grid_sampled(gid_list, input_type, stim_type, model_type, inputs, trial, 
									  fq_list, d_list, n_sample=3000, m=300):
    spikes_data = pd.DataFrame(columns=['fq', 'distance', 'p_value', 'v_length', 'v_angle', 
    	                                'n_spikes', 'v_length_std'])
    
    # Makes a grid to place the polar histograms in
    fig, axs = plt.subplots(len(fq_list), len(d_list), figsize=(20, 12), facecolor='w', edgecolor='k', subplot_kw=dict(polar=True))
    fig.subplots_adjust(hspace = 1.05, wspace=.0001)
    if len(fq_list) or len(d_list) > 1:
        axs = axs.ravel()
    
    # Generate a single concatenated table for all the cells in gid_list with the spike analysis columns
    table = generate_spike_analysis_table(gid_list, input_type, stim_type, model_type, inputs, trial)
    
    for i, (fq,d) in enumerate([[fq,d] for fq in fq_list for d in d_list]):
        all_spikes = []
        for index, row in table[(table['fq'] == fq) & (table['distance'] == d)].iterrows():
            if len(row['spike_phase_correct_window']) != 0:
                all_spikes = all_spikes + row['spike_phase_correct_window']

        if len(all_spikes) == 0:
            print 'There were no spikes at distance {} and frequency {}.'.format(d, fq)
            continue
            
        # Randomly sample n_sample spikes from all spikes
        spike_sample = random.sample(all_spikes, n_sample)

        v_angle = pycircstat.descriptive.mean(np.array(spike_sample))
        v_length = pycircstat.descriptive.vector_strength(np.array(spike_sample))
        p_value = pycircstat.rayleigh(np.array(spike_sample))[0]
        
        Y, X = np.histogram([x * 180./np.pi for x in spike_sample], bins=15)
        Xp =(X[1:] + X[:-1]) / 2
        Xp = Xp * np.pi / 180
        bars = axs[i].bar(Xp, Y, width=0.35)

        axs[i].set_ylim(0, m)
        axs[i].annotate('',xy=(v_angle, v_length * m), xytext=(v_angle,0), xycoords='data', arrowprops=dict(width=5, color='red'))
        axs[i].set_xlabel('Frequency: {}  Distance: {}\nNumber of Spikes: {}\nVector length: {:0.3f}\nVector angle: {:0.3f}\nP-value: {:0.3E}'.format(fq, d, len(spike_sample), v_length, v_angle * 180./np.pi, p_value))        
    
        spikes_data = spikes_data.append({'fq': fq, 
                                          'distance' : d, 
                                          'p_value' : p_value, 
                                          'v_length' : v_length, 
                                          'v_angle' : v_angle * 180./np.pi, 
                                          'n_spikes' : len(spike_sample),
                                          'v_length_std': pycircstat.descriptive.std(np.array(spike_sample))}, ignore_index=True)
    plt.tight_layout()    
    return axs, spikes_data