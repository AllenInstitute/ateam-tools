######################################################
# Authors: Fahimeh Baftizadeh, Taylor Connington
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
from neuroanalysis.miesnwb import MiesNwb
# import isee_engine.bionet.config as config
import table_plot_helper as tpl

"""
Many utils. Most are for building paths or file names or dealing with IO.
These functions solidify the contract between the data and the experimental/analysis and are used everywhere that data
is used--so best to be careful when changing the contract.
Ideally utils are as small as possible--a single unit of work--so that you can easily compose them
"""
################################################
#
#     Enumerated constant values
#
#################################################

# Enumerate stim types to use as source of truth
# Note: despite what it may look like elsewhere StimType.DC != 'dc' (!)
# Enumerate model types to use as source of truth

class Name(Enum):
    def __str__(self):
        return str(self.value)

class ElectrodeType(Name):
    EXTRA = 'extracellular'
    INTRA = 'intracellular'
    STIM = 'stimulus'

#################################################
#
#     Paths / filenames
#
#################################################

def concat_path(*args):
    """ Join paths together by parts, worry-free """
    root = args[0]
    is_abs_path = root[0] == '/'
    clean = [str(s).strip('/') for s in args]
    if is_abs_path:
        clean[0] = '/' + clean[0]
    return '/'.join(clean)

def get_dir_root(saved_data):
    """ For dealing with different netweork locations on Mac/Linux """
    network_root = '/allen'
    if os.path.isdir('/Volumes'):
        network_root = '/Volumes'
    if saved_data:
        network_root ='/local2'
    return network_root

def get_exp_id(conf):
    return conf['experiment info']['id']

def get_sampling_freq(conf):
    return conf['experiment info']['sampling freq']

def get_el_prop(conf, el_type):
    if el_type == ElectrodeType.EXTRA:
        return conf['extracellular electrode']
    if el_type == ElectrodeType.INTRA:
        return conf['intracellular electrode']
    if el_type == ElectrodeType.STIM:
        return conf['stimulus electrode']

def get_nwb_dir(conf):
    return conf['manifest']['$NWB_DIR']

def get_nwb_filename(exp_id, sampling_freq):
    return '{}_{}_compressed.nwb'.format(exp_id,sampling_freq)

def get_config_resolved_path(exp_id, sampling_freq, saved_data):
    network_root = get_dir_root(saved_data)
    return concat_path(network_root, 'aibs/mat/sooyl/result_tables',   str(exp_id) + "_" + str(sampling_freq) + '_resolved.json')


def get_config_path(exp_id, sampling_freq, saved_data): #Sampling frequency should be given in KHz
    """ Get dir containing runs for given params """
    network_root = get_dir_root(saved_data)
    keys = resolve_global_id(exp_id, sampling_freq)
    keys = keys + '.json'
    return concat_path(network_root, 'aibs/mat/sooyl/Stimulus_Item_Values', keys)

def resolve_global_id(exp_id, sampling_freq):
    parts = [exp_id , str(sampling_freq)]
    return "_".join(parts)

def get_table_dir(filename, saved_data):
    root_dir = get_dir_root(saved_data)
    return concat_path(root_dir,'aibs/mat/sooyl/result_tables/', filename)

def get_table_filename(exp_id, sampling_freq):
    return 'table_{}_{}.h5'.format(exp_id, sampling_freq)


def resolve_nwb_folder(sampling_freq):
    parts = [str(sampling_freq) + 'KHz_nwb_files' ]
    return "".join(parts)

def get_nwb_path_from_exp_id(exp_id, sampling_freq, saved_data):
    filename = get_nwb_filename(exp_id,sampling_freq)
    folder = resolve_nwb_folder(sampling_freq)
    root_dir = get_dir_root(saved_data)
    return concat_path(root_dir, 'aibs/mat/sooyl/',folder ,filename)

#################################################
#
#     Run aggregate data output file
#
#################################################

basic_cols = ['sweep_number','ex_el_id', 'stim_el_id', 'in_el_id', 'ex_el_distance(mu)', 'spike_tt']
ex_stimulus_cols = ['ex_amp(nA)', 'ex_frequency', 'ex_dur(ms)', 'ex_delay(ms)']
in_stimulus_cols = ['in_amp(pA)', 'in_dur(ms)', 'in_delay(ms)']
vi_phase_analysis_cols = ['vi_amp(mV)', 'vi_phase']
vext_phase_analysis_cols = ['vext_amp(mV)', 'vext_phase']
vm_phase_analysis_cols = ['vm_amp(mV)', 'vm_phase', 'avg_vm(mV)']
spike_phase_analysis_cols = ['spike_phase']

def resolve_additional_cols(include_ex_stimulus, include_in_stimulus, include_vext_phase_analysis,
                            include_vi_phase_analysis, include_spike_phase, include_vm_phase_analysis):
    """ Find all the additional columns for the table"""
    cols = basic_cols

    if include_ex_stimulus:
        cols= cols + ex_stimulus_cols
    if include_in_stimulus:
        cols= cols + in_stimulus_cols
    if include_vext_phase_analysis:
        cols = cols + vext_phase_analysis_cols
    if include_vi_phase_analysis:
        cols = cols + vi_phase_analysis_cols
    if include_spike_phase:
        cols = cols + spike_phase_analysis_cols
    if include_vm_phase_analysis:
        cols = cols + vm_phase_analysis_cols

    return cols

def resolve_run_id(sweep, exel, ex_amp, ex_freq, in_amp):
    stringified = map(str, [sweep, exel ,int(ex_amp), int(ex_freq), int(in_amp)])
    return '_'.join(stringified)



def build_df(additional_cols=[]):
    """ Wrapped df creation to give place to explicitly declare column types """
    df = pd.DataFrame(columns=(additional_cols))
    # pd doesn't do a great job of identifying ints
    # df['trial'] = df['trial'].astype(int)
    # df['electrode'] = df['electrode'].astype(int)

    return df

#################################################
#
#     Read / write
#
#################################################

def write_table_h5(fpath, df, attrs=None):
    """ dataframe to h5 """
    df.sort_values('sweep_number', inplace=True)
    with h5.File(fpath, 'w') as f5:
        f5.create_dataset('ids', data=map(str, df.index))
        for col in df.columns:
            if col not in ['spike_tt', 'spike_phase']: # cuz it will explode if it is type 'object'
                f5.create_dataset(col, data=df[col].astype(float))

        if set(['spike_tt', 'spike_phase']).issubset(df.columns):
            spike_threshold_t_grp = f5.create_group('spike_tt')
            spike_phase_grp = f5.create_group('spike_phase')

            for (rid, s_t_t) in df.spike_tt.iteritems():
                spike_threshold_t_grp.create_dataset(rid, maxshape=(None,), chunks=True, data=s_t_t)
        #
            for (rid, s_p) in df.spike_phase.iteritems():
                spike_phase_grp.create_dataset(rid, maxshape=(None,), chunks=True, data=s_p)

        if attrs is not None:
            for k,v in attrs.iteritems():
                # print k, v
                f5.attrs[k] = v



def read_table_h5(fpath):
    """ h5 to dataframe """

    with h5.File(fpath, 'r') as f5:

        extra_cols = []
        if 'has_ex_stimulus_data' in f5.attrs:
            extra_cols = extra_cols + ex_stimulus_cols

        if 'has_in_stimulus_data' in f5.attrs:
            extra_cols = extra_cols + in_stimulus_cols

        if 'include_vext_phase_analysis' in f5.attrs:
            extra_cols = extra_cols + vext_phase_analysis_cols

        if 'include_vi_phase_analysis' in f5.attrs:
            extra_cols = extra_cols + vi_phase_analysis_cols

        if 'include_vm_phase_analysis' in f5.attrs:
            extra_cols = extra_cols + vm_phase_analysis_cols

        if 'include_spike_phase' in f5.attrs:
            extra_cols = extra_cols + spike_phase_analysis_cols

        table = build_df(basic_cols + extra_cols)
        un_touched_data_cols = (basic_cols + extra_cols)
        ids = f5['ids'].value

        if set(['spike_tt', 'spike_phase']).issubset(un_touched_data_cols):
            spike_threshold_t_data = f5['spike_tt']
            spike_phase_data = f5['spike_phase']

        for i, rid in enumerate(ids):  # i is key for f5 dsets, rid is key for spikes & df
            data_cols = (basic_cols + extra_cols)

            if set(['spike_tt','spike_phase']).issubset(data_cols):
                spike_threshold_t_index = data_cols.index('spike_tt')
                spike_phase_index = data_cols.index('spike_phase')


            if set(['spike_tt', 'spike_phase']).issubset(data_cols):
                new_spike_threshold_t_index = data_cols.index('spike_tt')
                data_cols.pop(new_spike_threshold_t_index)
                new_spike_phase_index = data_cols.index('spike_phase')
                data_cols.pop(new_spike_phase_index)

            # print data_cols
            data = [f5[c][i] for c in data_cols]
            # place spikes in correct position
            if set(['spike_tt','spike_phase']).issubset(un_touched_data_cols):
                data.insert(spike_threshold_t_index, spike_threshold_t_data[rid].value)
                data.insert(spike_phase_index, spike_phase_data[rid].value)

            table.loc[rid] = data
        table['spike_tt_A'] = np.NAN
        table['spike_phase_A'] = np.NAN
        table['spike_tt_A'] = table['spike_tt_A'].astype(object)
        table['spike_phase_A'] = table['spike_phase_A'].astype(object)

        for index, row in table.iterrows():
            table['spike_tt_A'][index] = tpl.filter_list(row['in_amp(pA)'], row['ex_delay(ms)'], row['ex_dur(ms)'],
                                                     row['spike_tt'], row['spike_tt'])
            table['spike_phase_A'][index] = tpl.filter_list(row['in_amp(pA)'], row['ex_delay(ms)'], row['ex_dur(ms)'],
                                                        row['spike_tt'], row['spike_phase'])

        table['num_spikes'] = table.apply(lambda row: len(row['spike_tt']), axis=1)
        table['num_spikes_A'] = table.apply(lambda row: len(row['spike_tt_A']), axis=1)
        table['spike_phase_A_corrected'] = table['spike_phase_A'].apply(tpl.phase_correction)


        # table['num_spikes'] = table.apply(lambda row: len(row['spike_tt']), axis=1)
        # table['spike_tt_A'] = table.apply(tpl.filter_list, varname1="spike_tt", varname2="spike_tt", axis=1)
        # table['spike_phase_A'] = table.apply(tpl.filter_list, varname1="spike_tt", varname2="spike_phase", axis=1)
        # table['num_spikes_A'] = table.apply(lambda row: len(row['spike_tt_A']), axis=1)
        # table['spike_phase_A_corrected'] = table.apply(tpl.phase_correction, axis=1)
    return table



def read_table_from_exp_id(exp_id, sampling_freq, saved_data):
    filename = get_table_filename( exp_id, sampling_freq)
    dir = get_table_dir(filename, saved_data)
    table = read_table_h5(dir)
    return table

def read_trace_from_nwb(exp_id, sampling_freq,  sweep, el_id, saved_data):
    nwb_path = get_nwb_path_from_exp_id(exp_id, sampling_freq, saved_data)
    resolved_config_path = get_config_resolved_path(exp_id, sampling_freq, saved_data)
    nwb = MiesNwb(nwb_path)
    # conf = config.build(resolved_config_path)
    #stimulus_description = nwb.contents[sweep][el_id].stimulus.description
    #BB = conf['stimulus description'][str(stimulus_description)]
    v = nwb.contents[sweep][el_id]['primary'].data * 1000
    return  v













# def format_sweep(sweep_number):
#     """ Formatting for el number for cleaner organization. Idempotent """
#     return str(sweep_number).zfill(5)
#
# 
# def get_param_dir(saved_data, *args):
#     network_root = get_dir_root(saved_data)
#     return concat_path(network_root, 'aibs/mat/sooyl/parameter_sheets/', *args)
#
# # def get_nwb_dir(sampling_freq, saved_data, *args): #Sampling frequency should be given in KHz
# #     """ Get dir containing runs for given params """
# #     network_root = get_dir_root(saved_data)
# #     nwb_folder = resolve_nwb_folder(sampling_freq)
# #     return concat_path(network_root, 'aibs/mat/sooyl/', nwb_folder, *args)
#
# def get_param_filename(exp_id,sampling_freq):
#     return '{}_{}_parameters.xlsx'.format(exp_id, sampling_freq)
#
#
#
# def get_nwb_file_path(exp_id, sampling_freq, saved_data):
#     parts = [get_nwb_dir(sampling_freq, saved_data) +'/' + get_nwb_filename(exp_id,sampling_freq)]
#     return "".join(parts)
#
# def get_param_file_path(exp_id, sampling_freq, saved_data):
#     parts = [get_param_dir(saved_data) +'/' + get_param_filename(exp_id,sampling_freq)]
#     return "".join(parts)
#
#
#
# def get_electrode_id(which_electrode, exp_id, sampling_freq, saved_data):
#     param_df = read_param_file(exp_id, sampling_freq, saved_data)
#     electrode_id = [x.encode('ascii','ignore') for x  in param_df[which_electrode].unique().tolist() if x is not np.nan]
#     if len(electrode_id) > 1:
#         print "ERROR, in the electrode id column, there should be only one label"
#     else:
#         return electrode_id[0]
#
# def get_sweep_column_name(sweep_number, which_electrode, exp_id, sampling_freq, saved_data):
#     electrode_id = get_electrode_id(which_electrode, exp_id, sampling_freq, saved_data)
#     parts = ['data' , format_sweep(sweep_number) , electrode_id]
#     sweep_col_name =  '_'.join(parts)
#     parts = [ 'acquisition/timeseries',sweep_col_name, 'data']
#     return '/'.join(parts)
#
#
#
# def read_nwb_file(exp_id, sampling_freq, saved_data):
#     path = get_nwb_file_path(exp_id, sampling_freq, saved_data)
#     cvh5 = h5.File(path)
#     return cvh5
#
# def read_param_file(exp_id, sampling_freq, saved_data):
#     path = get_param_file_path(exp_id, sampling_freq, saved_data)
#     return pd.read_excel(path)
#
# def read_timeseries(sweep_number, which_electrode, exp_id, sampling_freq, saved_data):
#     cvh5 = read_nwb_file(exp_id, sampling_freq, saved_data)
#     sweep_column_name = get_sweep_column_name(sweep_number, which_electrode, exp_id, sampling_freq, saved_data)
#     return cvh5[sweep_column_name].value
