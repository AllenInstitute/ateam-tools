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

"""
Many utils. Most are for building paths or file names or dealing with IO. 
These functions solidify the contract between the data and the simulation/analysis and are used everywhere that data 
is used--so best to be careful when changing the contract. 
Ideally utils are as small as possible--a single unit of work--so that you can easily compose them.
Names:
get...dir  -> returns path to a directory. last param *args allows additional paths to be added, with the exception of
            specific dc/sin/etc output dirs, since these are basically endpoints
get...path -> returns path to a file
resolve... -> the function is essentially a definition for whatever it outputs
"""


#################################################
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

class StimType(Name):
    DC = 'dc'
    SIN = 'sin'
    SIN_DC = 'sin_dc'
    SIN_DC_LGN = 'sin_dc_lgn'
    SIN_DC_POISSON = 'sin_dc_poisson'

class InputType(Name):
    EXTRASTIM = 'extrastim'
    EXTRASTIM_INTRASTIM = 'extrastim_intrastim'
    EXTRASTIM_INTRASTIM_SYN = 'extrastim_intrastim_syn'
    INTRASTIM = "intrastim"

class ModelType(Name):
    PERISOMATIC = 'perisomatic'
    ACTIVE      = 'all_active'
    PASSIVE     = 'passive'
    FAHIMEH    = 'fahimeh'

#################################################
#
#     Manage arguments
#
#################################################

def get_input_args(cell_id, input_type, stim_type, model_type, trial=0, **kwargs):
    if 'el' and 'amp' not in kwargs.keys():
        print ("el or amp are not included in the arguments")

    if input_type == "extrastim_intrastim":
        if stim_type != "sin_dc":
            raise ValueError("stim_type should match input_type")
        if 'freq' and 'ic_amp' not in kwargs.keys():
            raise ValueError("You are required to include freq and ic_amp in the arguments for this input_type")

    if input_type == "extrastim":
        if stim_type != "sin":
            raise ValueError("stim_type should match input_type")
        if 'freq' not in kwargs.keys():
            raise ValueError("You are required to include freq in the arguments for this input_type")
        if 'ic_amp' in kwargs.keys():
            raise ValueError("Ic_amp should not be in the arguments for this input_type")

    el = kwargs['el'] if 'el' in kwargs else None
    amp = kwargs['amp'] if 'amp' in kwargs else None
    freq = kwargs['freq'] if 'freq' in kwargs else None
    ic_amp = kwargs['ic_amp'] if 'ic_amp' in kwargs else None

    return el, amp, freq, ic_amp


#################################################
#
#     Paths / filenames
#
#################################################

def format_el(el):
    """ Formatting for el number for cleaner organization. Idempotent """
    return str(el).zfill(4)

def format_amp(amp):
    """ Formatting of amp for easier naming. Idempotent """
    return amp if type(amp) == str else "{0:.0f}".format( math.fabs(amp * 1000000.)) # takes in mA and Returns in nA

def format_icamp(ic_amp):
    """ Formatting of amp for easier naming. Idempotent """
    return ic_amp if type(ic_amp) == str else "{0:.0f}".format( ic_amp * 1000.) # takes in nA and Return in PA

def format_freq(freq):
    """ Formatting of freq for naming. Idempotent. """
    return freq if type(freq) == str else "{0:.0f}".format(freq)

def get_dir_root(saved_data):
    """ For dealing with different netweork locations on Mac/Linux """
    network_root = '/allen'
    if os.path.isdir('/Volumes'):
        network_root = '/Volumes'
    if saved_data:
        network_root ='/local2'
    return network_root

def concat_path(*args):
    """ Join paths together by parts, worry-free """
    root = args[0]
    is_abs_path = root[0] == '/'
    clean = [str(s).strip('/') for s in args]
    if is_abs_path:
        clean[0] = '/' + clean[0]
    return '/'.join(clean)

def get_reobase_dir(saved_data, *args):
    network_root = get_dir_root(saved_data)
    return concat_path(network_root, 'aibs/mat/Fahimehb/Data_cube/reobase', *args)
    # return concat_path('/allen/aibs/mat/maddie/Bionet_example/', *args)

def get_output_dir(input_type, stim_type, model_type, cell_gid, saved_data, *args):
    """ Get dir containing runs for given params """
    reobase_dir = get_reobase_dir(saved_data)
    return concat_path(reobase_dir, 'Run_folder/outputs/', input_type, stim_type, model_type, cell_gid, *args)

def get_intrastim_output_dir(cell_gid, input_type, stim_type, model_type, ic_amp, trial, saved_data):
    root_dir = get_output_dir(input_type , stim_type, model_type, cell_gid, saved_data)
    out_dir = get_intrastim_dir_name(ic_amp, trial)
    return concat_path(root_dir, out_dir)

def get_electrode_path(electrodes_dir, gid, el):
    return concat_path(electrodes_dir, str(gid) + '_' + format_el(el) + '.csv')

def get_config_resolved_path(stim_type, out_folder, el, amp, freq= None, ic_amp= None):
    if stim_type == StimType.DC:
        key = resolve_dc_key(el, amp)
    if stim_type == StimType.SIN:
        key = resolve_sin_key(el, amp, freq)
    if stim_type == StimType.SIN_DC:
        key = resolve_sin_dc_key(el,amp, freq, ic_amp)

    return concat_path(out_folder, 'config_' + key + '_resolved.json')

### reolve the stim_type

def resolve_key(el, amp, freq=None, ic_amp=None):
    if all(v is None for v in [freq, ic_amp]):
        resolve_dc_key(el,amp)
    elif freq is not None and ic_amp is None:
        resolve_sin_key(el, amp, freq)
    else:
        resolve_sin_dc_key(el, amp, freq, ic_amp)

def get_dir_name(el, amp, freq=None, ic_amp=None, trial =0):
    if all(v is None for v in [freq, ic_amp]):
        return get_dc_dir_name(el, amp, trial)
    elif freq is not None and ic_amp is None:
        return get_sin_dir_name(el, amp, freq, trial)
    else:
        return get_sin_dc_dir_name(el, amp, freq, ic_amp, trial)


# def get_outout_dir(cell_gid, input_type, stim_type, model_type, el, amp, trial= 0 , freq=None, ic_amp = None):
#     if stim_type == "dc":
#         get_dc_output_dir(cell_gid, input_type, stim_type, model_type, el, amp, trial)
#     if stim_type == "sin":
#         get_sin_output_dir(cell_gid, input_type, stim_type, model_type, el, amp, freq, trial)
#     if stim_type == "sin_dc":
#         get_sin_dc_output_dir(cell_gid, input_type, stim_type, model_type, el, amp, freq, ic_amp, trial)

### INTRASTIM_DC ###

def resolve_intrastim_dc_key(icamp):
    parts = ['icamp' + format_icamp(icamp)]
    return "_".join(parts)

def get_intrastim_dir_name(icamp, trial):
    return '_'.join([resolve_intrastim_dc_key(icamp), 'tr' + str(trial)])


### SIN_DC ###

def resolve_sin_dc_key(el, amp, freq, ic_amp):
    """ file/folder name code """
    parts = ['el' + format_el(el), 'amp' + format_amp(amp), 'freq' + format_freq(freq), 'icamp' + format_icamp(ic_amp)]
    return '_'.join(parts)

def get_sin_dc_dir_name(el, amp, freq, ic_amp, trial):
    return '_'.join([resolve_sin_dc_key(el, amp, freq, ic_amp), 'tr' + str(trial)])

# def get_sin_dc_output_dir(cell_gid, input_type, stim_type, model_type, el, amp, freq, ic_amp, trial=0):
#     root_dir = get_output_dir(input_type , stim_type, model_type, cell_gid)
#     out_dir = get_sin_dc_dir_name(el, amp, freq, ic_amp, trial)
#     return concat_path(root_dir, out_dir)


### DC ###

def resolve_dc_key(el, amp):
    """ file/folder name code """
    parts = ['el' + format_el(el), 'amp' + format_amp(amp)]
    return '_'.join(parts)

def get_dc_dir_name(el, amp, trial):
    return '_'.join([resolve_dc_key(el, amp), 'tr' + str(trial)])

def get_dc_output_dir(cell_gid, input_type, stim_type, model_type, el, amp, trial):
    root_dir = get_output_dir(input_type, stim_type, model_type, cell_gid)
    out_dir = get_dc_dir_name(el, amp, trial)
    return concat_path(root_dir, out_dir)

### SIN ###

def resolve_sin_key(el, amp, freq):
    """ file/folder name code """
    parts = ['el' + format_el(el), 'amp' + format_amp(amp), 'freq' + format_freq(freq)]
    return '_'.join(parts)

def get_sin_dir_name(el, amp, freq, trial):
    return '_'.join([resolve_sin_key(el, amp, freq), 'tr' + str(trial)])

def get_sin_output_dir(cell_gid, input_type, stim_type, model_type, el, amp, freq, trial=0):
    root_dir = get_output_dir(input_type, stim_type, model_type, cell_gid)
    out_dir = get_sin_dir_name(el, amp, freq, trial)
    return concat_path(root_dir, out_dir)

### Tables ###

def get_table_dir(input_type, stim_type, model_type, *args):
    # print concat_path('/allen/aibs/mat/Fahimehb/Data_cube/reobase/Run_folder/result_tables/', input_type, stim_type, model_type, *args)
    return concat_path('/allen/aibs/mat/Fahimehb/Data_cube/reobase/Run_folder/result_tables/', input_type, stim_type, model_type, *args)
    # return concat_path('/allen/aibs/mat/maddie/Bionet_example/Run_folder/result_tables/', input_type, stim_type, model_type, *args)

def get_table_filename(cell_gid, amp,trial):
    return 'table_{}_amp{}_tr{}.h5'.format(cell_gid, format_amp(amp), trial)

def get_control_table_filename(cell_gid, amp,trial):
    return 'control_table_{}_amp{}_tr{}.h5'.format(cell_gid, format_amp(amp), trial)

### Modulation_tables ###
def get_analysis_dir(input_type, stim_type, model_type, *args):
    # print concat_path('/allen/aibs/mat/Fahimehb/Data_cube/reobase/Run_folder/result_tables/', input_type, stim_type, model_type, *args)
    return concat_path('/allen/aibs/mat/Fahimehb/Data_cube/reobase/Run_folder/Analysis/result_modulation_tables/', input_type, stim_type, model_type, *args)
    # return concat_path('/allen/aibs/mat/maddie/Bionet_example/Run_folder/result_modulation_tables/', input_type, stim_type, model_type, *args)

def get_modulation_table_filename(cell_gid,trial):
    return '{}_vm_amp_modulation_tr{}.csv'.format(cell_gid, trial)

## VMDs ###

def get_vmd_dir(input_type, stim_type, model_type, *args):
    return get_reobase_dir(None, 'Run_folder/result_vmd/', input_type, stim_type, model_type, *args)

def get_vmd_filename(cell_gid, amp, vmdtype, trial):
    return '{}_amp{}_{}_tr{}.pdb'.format(cell_gid, format_amp(amp), vmdtype, trial)

#################################################
#
#     IO / data
#
#################################################

def get_spikes_file(output):
    return h5.File(concat_path(output, 'spikes.h5'), 'r')


def get_cv_files(output, cells=None):
    cv_dir = concat_path(output, 'cellvars')

    if cells is not None:
        files = map(lambda c: h5.File(concat_path(cv_dir, str(c)) + '.h5', 'r'), cells)
    else:
        files = [h5.File(f) for f in glob.glob(cv_dir + "/*.h5")]

    return files

def get_cv_files_f_i(output, current, cells=None):
    # Just to get cv files for f-I plots
    # TODO: improve plot_f_i function or get_cv_files function so can use get_cv_files directly
    cv_dir = concat_path(output, 'el0000_amp{}_tr0'.format(current), 'cellvars')

    if cells is not None:
        files = map(lambda c: h5.File(concat_path(cv_dir, str(c)) + '.h5', 'r'), cells)
    else:
        files = [h5.File(f) for f in glob.glob(cv_dir + "/*.h5")]

    return files

def get_json_from_file(path):
    with open(path, 'r') as f:
        return json.load(f)

def get_spikes_data(out_path):
    cellvars = get_cv_files(out_path, cells=[0])[0]
    return cellvars['spikes'].value

def get_spikes_data2(out_path):
    ## UNTESTED! hoping this will be faster than using cv files.......
    ## TODO finish
    spikes = pd.read_csv(concat_path(out_path,'spikes.txt'), ' ')
    cellvars = get_cv_files(out_path, cells=[0])[0]
    return cellvars['spikes'].value

def get_electrode_xyz(electrode_pos_path):  # Ideally you would use the method bionet uses
    # mesh files are unnecessary for this study
    electrode_pos_df = pd.read_csv(electrode_pos_path, sep=' ')
    return electrode_pos_df.as_matrix(columns=['pos_x', 'pos_y', 'pos_z'])


def get_cell_xyz(cell_file):  # Ideally you would use the method bionet uses
    cell_props_df = pd.read_csv(cell_file, sep=' ')
    return cell_props_df.as_matrix(columns=['x_soma', 'y_soma', 'z_soma'])


#################################################
#
#     Run aggregate data output file
#
#################################################

dc_cols = ['trial', 'electrode', 'x', 'y', 'z', 'distance', 'amp', 'spikes']
vm_cols = ['vm_stim', 'delta_vm']
iclamp_cols = ['ic_amp']
sin_cols=['fq']
vm_phase_analysis_cols = ['vm_amp', 'vm_phase']
vext_phase_analysis_cols = ['vext_amp', 'vext_phase']
vi_phase_analysis_cols = ['vi_amp', 'vi_phase']
spike_phase_analysis_cols = ['spike_threshold_t', 'spike_phase']
sta_cols = ['sta']

def resolve_additional_cols(include_delta_vm, include_sin, include_iclamp, include_vm_phase_analysis,
                            include_vext_phase_analysis, include_vi_phase_analysis, include_spike_phase, include_sta):
    """ Find all the additional columns for the table"""
    cols = []

    if include_delta_vm:
        cols= cols + vm_cols
    if include_sin:
        cols= cols + sin_cols
    if include_iclamp:
        cols= cols + iclamp_cols
    if include_vm_phase_analysis:
        cols= cols + vm_phase_analysis_cols
    if include_vext_phase_analysis:
        cols = cols + vext_phase_analysis_cols
    if include_vi_phase_analysis:
        cols = cols + vi_phase_analysis_cols
    if include_spike_phase:
        cols = cols + spike_phase_analysis_cols
    if include_sta:
        cols = cols + sta_cols

    return cols

def resolve_output_aggregates(include_delta_vm, include_sin, include_iclamp, amp, trial):
    """ Find all the outputs for the table"""
    if  not include_sin and not include_iclamp:
        return get_dc_dir_name(9999, amp, trial).replace('9999', '*')
    if include_sin and not include_iclamp:
        return get_dir_name(el=9999, amp=amp, freq=999, ic_amp=None, trial=trial).replace('9999', '*').replace('999', '*')
    if include_iclamp:
        return get_dir_name(el=9999, amp=amp, freq=999, ic_amp=0.099, trial=trial).replace('9999', '*').replace('999', '*').replace('99', '*')


# def resolve_dc_run_id(gid, electrode, amp):
#     """ Unique run id (per stim type) """
#     stringified = map(str, [gid, electrode, format_amp(amp)])
#     return '_'.join(stringified)  # using current in micro amps
#
# def resolve_sin_dc_run_id(gid, electrode, amp, freq, ic_amp):
#     stringified = map(str, [gid, electrode, format_amp(amp), freq, ic_amp])
#     return '_'.join(stringified)
#
# def resolve_sin_run_id(gid, electrode, amp, freq):
#     stringified = map(str, [gid, electrode, format_amp(amp), freq])
#     return '_'.join(stringified)

def resolve_run_id(gid, electrode, amp, freq=None, ic_amp=None):
    if freq is not None and ic_amp is not None:
        stringified = map(str, [gid, electrode, format_amp(amp), freq, format_icamp(ic_amp)])
    if freq is not None and ic_amp is None:
        stringified = map(str, [gid, electrode, format_amp(amp), freq])
    if freq is None and ic_amp is None:
        stringified = map(str, [gid, electrode, format_amp(amp)])
    return '_'.join(stringified)


def build_df(additional_cols=[]):
    """ Wrapped df creation to give place to explicitly declare column types """
    df = pd.DataFrame(columns=(additional_cols))
    # pd doesn't do a great job of identifying ints
    df['trial'] = df['trial'].astype(int)
    df['electrode'] = df['electrode'].astype(int)

    return df


def read_table_h5(fpath):
    """ h5 to dataframe """

    with h5.File(fpath, 'r') as f5:

        extra_cols = []
        if 'has_vm_data' in f5.attrs and f5.attrs['has_vm_data']:
            extra_cols = extra_cols + vm_cols

        if 'has_sin' in f5.attrs and f5.attrs['has_sin']:
            extra_cols = extra_cols + sin_cols

        if 'has_iclamp' in f5.attrs and f5.attrs['has_iclamp']:
            extra_cols = extra_cols + iclamp_cols

        if 'has_vm_phase_analysis' in f5.attrs and f5.attrs['has_vm_phase_analysis']:
            extra_cols = extra_cols + vm_phase_analysis_cols

        if 'has_vext_phase_analysis' in f5.attrs and f5.attrs['has_vext_phase_analysis']:
            extra_cols = extra_cols + vext_phase_analysis_cols

        if 'has_vi_phase_analysis' in f5.attrs and f5.attrs['has_vi_phase_analysis']:
            extra_cols = extra_cols + vi_phase_analysis_cols

        if 'has_spike_phase_analysis' in f5.attrs and f5.attrs['has_spike_phase_analysis']:
            extra_cols = extra_cols + spike_phase_analysis_cols

        # if 'has_sta_analysis' in f5.attrs and f5.attrs['has_sta_analysis']:
        #     extra_cols = extra_cols + sta_cols

        table = build_df(dc_cols + extra_cols)
        un_touched_data_cols = (dc_cols + extra_cols)
        ids = f5['ids'].value

        if set(['spike_threshold_t', 'spike_phase','spikes']).issubset(un_touched_data_cols):
            spike_threshold_t_data = f5['spike_threshold_t']
            spike_phase_data = f5['spike_phase']
            spikes_data = f5['spikes']

        for i, rid in enumerate(ids):  # i is key for f5 dsets, rid is key for spikes & df
            data_cols = (dc_cols + extra_cols)
            if set(['spike_threshold_t','spike_phase', 'spikes']).issubset(data_cols):
                spike_threshold_t_index = data_cols.index('spike_threshold_t')
                spike_phase_index = data_cols.index('spike_phase')
                spikes_index = data_cols.index('spikes')

            if set(['spike_threshold_t','spike_phase', 'spikes']).issubset(data_cols):
                new_spike_threshold_t_index = data_cols.index('spike_threshold_t')
                data_cols.pop(new_spike_threshold_t_index)
                new_spike_phase_index = data_cols.index('spike_phase')
                data_cols.pop(new_spike_phase_index)
                new_spikes_index = data_cols.index('spikes')
                data_cols.pop(new_spikes_index)

            # print   f5['spike_threshold_t'][i]
            # print data_cols
            data = [f5[c][i] for c in data_cols]
        #         # place spikes in correct position
            if set(['spike_threshold_t', 'spike_phase','spikes']).issubset(un_touched_data_cols):
                    data.insert(spike_threshold_t_index, spike_threshold_t_data[rid].value)
                    data.insert(spike_phase_index, spike_phase_data[rid].value)
                    data.insert(spikes_index, spikes_data[rid].value)

            table.loc[rid] = data
        return table
        # table = build_df(extra_cols)
        # un_touched_data_cols = (dc_cols + extra_cols)
        # print dc_cols + extra_cols
        # ids = f5['ids'].value
        # spike_data = f5['spikes']
        #
        # if set(['spike_threshold_t', 'spike_phase', 'sta']).issubset(un_touched_data_cols):
        #     spike_threshold_t_data = f5['spike_threshold_t']
        #     spike_phase_data = f5['spike_phase']
        #     sta_data = f5['sta']
        #
        # for i, rid in enumerate(ids):  # i is key for f5 dsets, rid is key for spikes & df
        #     data_cols = (dc_cols + extra_cols)
        #     spike_index = data_cols.index('spikes')
        #
        #     if set(['spike_threshold_t','spike_phase', 'sta']).issubset(data_cols):
        #         spike_threshold_t_index = data_cols.index('spike_threshold_t')
        #         spike_phase_index = data_cols.index('spike_phase')
        #         sta_index = data_cols.index('sta')
        #
        #     data_cols.pop(spike_index)
        #     if set(['spike_threshold_t', 'spike_phase', 'sta']).issubset(data_cols):
        #         new_spike_threshold_t_index = data_cols.index('spike_threshold_t')
        #         data_cols.pop(new_spike_threshold_t_index)
        #         new_spike_phase_index = data_cols.index('spike_phase')
        #         data_cols.pop(new_spike_phase_index)
        #         new_sta_index = data_cols.index('sta')
        #         data_cols.pop(new_sta_index)
        #
        #
        #     data = [f5[c][i] for c in data_cols]
        #     # place spikes in correct position
        #     data.insert(spike_index, spike_data[rid].value)
        #     if set(['spike_threshold_t','spike_phase', 'sta']).issubset(un_touched_data_cols):
        #         data.insert(spike_threshold_t_index, spike_threshold_t_data[rid].value)
        #         data.insert(spike_phase_index, spike_phase_data[rid].value)
        #         data.insert(sta_index, sta_data[rid].value)
        #     table.loc[rid] = data

    # return table


def write_table_h5(fpath, df, attrs=None):
    """ dataframe to h5 """
    df.sort_values('electrode', inplace=True)
    with h5.File(fpath, 'w') as f5:
        f5.create_dataset('ids', data=map(str, df.index))
        for col in df.columns:
            if col not in ['spikes', 'spike_threshold_t', 'spike_phase', 'sta']: # cuz it will explode if it is type 'object'
                f5.create_dataset(col, data=df[col])

        spike_grp = f5.create_group('spikes')
        for (rid, spike_times) in df.spikes.iteritems():
            spike_grp.create_dataset(rid, maxshape=(None,), chunks=True, data=spike_times)


        if set(['spike_threshold_t', 'spike_phase', 'sta']).issubset(df.columns):
            spike_threshold_t_grp = f5.create_group('spike_threshold_t')
            spike_phase_grp = f5.create_group('spike_phase')
            sta_grp = f5.create_group('sta')

            for (rid, s_t_t) in df.spike_threshold_t.iteritems():
                spike_threshold_t_grp.create_dataset(rid, maxshape=(None,), chunks=True, data=s_t_t)
        #
            for (rid, s_p) in df.spike_phase.iteritems():
                spike_phase_grp.create_dataset(rid, maxshape=(None,), chunks=True, data=s_p)

            for (rid, sta_p) in df.sta.iteritems():
                sta_grp.create_dataset(rid, maxshape=(None,), chunks=True, data=sta_p)

        if attrs is not None:
            for k,v in attrs.iteritems():
                # print k, v
                f5.attrs[k] = v


def read_cell_tables(cell_gid, amp_range, input_type, stim_type, model_type, trial,
                     control = False, data_dir=None):
    """ Read h5 files for a set of amplitudes """
    # print "Fetching data..."
    data_dir = get_table_dir(input_type, stim_type, model_type) if data_dir is None else data_dir
    if control:
        paths = [concat_path(data_dir, get_control_table_filename(cell_gid, a, trial)) for a in amp_range]
    else:
        paths = [concat_path(data_dir, get_table_filename(cell_gid, a, trial)) for a in amp_range]

    t = build_df(dc_cols)  # do this for code analysis
    t = t.append([read_table_h5(p) for p in paths])
    t['num_spikes'] = t.apply(lambda row: len(row['spikes']), axis=1)
    # t['num_true_spikes'] = np.where(t['num_spikes'] == 1, 0, t['num_spikes'])

    # print "Done"
    return t


def read_cell_rows(cell_gid, els, amps, stim_type, trial,
                   data_dir=None):
    """
    Read h5 files for a set of amps and electrodes -- faster when considering a small subset of electrodes
    Doesn't use read_table_h5--it actually interfaces with h5 file.
    """
    data_dir = get_table_dir(stim_type, ModelType.PERISOMATIC) if data_dir is None else data_dir
    els = [els] if type(els) is not list else els
    amps = [amps] if type(amps) is not list else amps # non-formatted
    data_cols = [x for x in dc_cols if x != 'spikes']
    table = build_df(dc_cols)

    for amp in amps:
        fpath = concat_path(data_dir, get_table_filename(cell_gid, amp, trial))

        with h5.File(fpath, 'r') as f5:
            ids = f5['ids'].value
            spike_data = f5['spikes']

            for el in els:
                rid = resolve_run_id(cell_gid,el,amp)
                i = np.argwhere(ids == rid)[0][0] # i is key for f5 dsets, rid is key for spikes & df
                data = [f5[c][i] for c in data_cols]
                table.loc[rid] = data + [(spike_data[rid].value)]

    table['num_spikes'] = table.apply(lambda row: len(row['spikes']), axis=1)

    return table

def get_index_close_els(cell_gid, input_type, stim_type, model_type, trial, saved_data):
    """Get list of bad electrodes for which the simulation was not finished"""
    output_dir = get_output_dir(input_type, stim_type, model_type, cell_gid, saved_data)
    # print output_dir
    index_list=[]
    for filename in os.listdir(output_dir):
        if ("tr" + str(trial)) in filename: 
            log_file= concat_path(output_dir, filename, "/log.txt")
                #if 'External electrode is too close' in open(log_file).read():
            if not 'Simulation Duration' in open(log_file).read():
                el = int([x for x in filename.split('_') if x.startswith('el')][0][2:])
                ic_amp = int([x for x in filename.split('_') if x.startswith('icamp')][0][5:]) * 0.001 if "icamp" in filename else None
                fq = int([x for x in filename.split('_') if x.startswith('freq')][0][4:]) if "freq" in filename else None
                amp = int([x for x in filename.split('_') if x.startswith('amp')][0][3:]) * 0.000001
                index_list.append(resolve_run_id(gid=cell_gid, electrode=el, amp=amp, freq=fq, ic_amp=ic_amp))
    return index_list

#############################################

## VMD ##

#############################################

def df_to_pdb(output, df, x_col, y_col, z_col, beta_col):
    with open(output, 'w') as f:
        df.index = np.arange(1,len(df)+1)
        for index, row in df.iterrows():
            f.write('%s %6s %s %11.3f%8.3f%8.3f %s %5.2f\n' % ("ATOM",index,
                                                                 " O   PRO A   1",
                                                                 row[x_col],
                                                                 row[y_col],
                                                                 row[z_col],
                                                                 " 1.00",
                                                                 row[beta_col]))


def get_swc_filename(cell_models_filepath, cell_gid):
    df = pd.read_csv(cell_models_filepath, sep = " ", header = "infer")
    index = np.where(df["model_id"] == cell_gid)[0][0]
    return df["morphology"][index]


def read_swc_file(swc_filename):
    swc_list = pd.read_csv(swc_filename, sep=" ", comment="#",
                           names=["id", "type", "x", "y", "z", "r", "pid"])
    soma_pos = swc_list[swc_list["type"] == 1][["x", "y", "z"]]
    swc_list = swc_list[["x", "y", "z", "r"]]
    return swc_list, soma_pos


def move_swc_to_origin(swc_list, soma_pos):
    x = swc_list["x"] - float(soma_pos["x"])
    y = swc_list["y"] - float(soma_pos["y"])
    z = swc_list["z"] - float(soma_pos["z"])
    r = swc_list["r"]
    swc_movedto_origin = pd.DataFrame({"x": x, "y": y, "z": z, "r": r})

    return swc_movedto_origin


def table_to_pdb(gids, inputs, input_type, stim_type, model_type, trial):
    for gid in gids:
        for input in inputs:
            t = read_cell_tables(gid, [input], input_type, stim_type, model_type, trial)
            grouped_df = t.groupby('amp')
            gb = grouped_df.groups
            if (all(t["num_spikes"] < 2) == True):
                vmd_filename = get_vmd_filename(gid, input, "sub", trial)
                vmd_result_dir = get_vmd_dir(input_type, stim_type, model_type, vmd_filename)
                for key, values in gb.iteritems():
                    df_to_pdb(vmd_result_dir, t.loc[values], "x", "y", "z", "delta_vm")
            else:
                vmd_filename = get_vmd_filename(gid, input, "supra", trial)
                vmd_result_dir = get_vmd_dir(input_type, stim_type, model_type, vmd_filename)
                for key, values in gb.iteritems():
                    df_to_pdb(vmd_result_dir, t.loc[values], "x", "y", "z", "num_spikes")


def swc_to_pdb(gids, input_type, stim_type, model_type, cell_models_file):
    for gid in gids:
        path = get_reobase_dir(None, "Run_folder", "network", cell_models_file)
        swc_filename = get_swc_filename(path, gid)
        swc_filepath = get_reobase_dir(None, "cell_models_layer4", "biophysical", "morphology",
                                      swc_filename)
        swc_list, soma_pos = read_swc_file(swc_filepath)
        swc_moved_to_origin = move_swc_to_origin(swc_list, soma_pos)
        swcpdb_filename = "swc" + str(gid) + ".pdb"

        swcpdb_filedir = get_vmd_dir(input_type, stim_type, model_type, swcpdb_filename)
        df_to_pdb(swcpdb_filedir, swc_moved_to_origin, "x", "y", "z", "r" )


def csv_to_pdb(csv_path, pdb_filedir):
    csv_values = pd.read_csv(csv_path, sep =" ", comment="T", names=["x", "y", "z"])
    empty = pd.DataFrame([0.0])
    df_to_pdb(pdb_filedir, csv_values, "x", "y", "z", empty[0])


# def make_dir_if_new(path):
#     try: 
#         os.makedirs(path)
#     except OSError:
#         if not os.path.isdir(path):
#             raise