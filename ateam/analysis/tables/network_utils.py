######################################################
# Authors: Fahimeh Baftizadeh
# Date created: 4/1/2018
######################################################


import numpy as np
import pandas as pd
import os, sys
from enum import Enum
import math
import h5py as h5
import glob
import json

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

class InputType(Name):
    EXTRASTIM_LGN = 'extrastim_lgn'
    EXTRASTIM = 'extrastim'

class ModelType(Name):
    PERISOMATIC = 'perisomatic'
    ACTIVE      = 'all_active'
    PERISOMATIC_LIF     = 'perisomatic_lif'
    ACTIVE_LIF    = 'all_active_lif'

#################################################
#
#     Paths / filenames
#
#################################################

def format_amp(amp):
    """ Formatting of amp for easier naming. Idempotent """
    return amp if type(amp) == str else "{0:.0f}".format( math.fabs(amp * 1000000.)) # takes in mA and Returns in nA

def format_freq(freq):
    """ Formatting of freq for naming. Idempotent. """
    return freq if type(freq) == str else "{0:.0f}".format(freq)

def get_dir_root(saved_data):
    """ For dealing with different netweork locations on Mac/Linux """
    network_root = '/allen/aibs'
    if os.path.isdir('/Volumes'):
        network_root = '/Volumes'
    if saved_data:
        network_root ='/local2/aibs'
    return network_root

def concat_path(*args):
    """ Join paths together by parts, worry-free """
    root = args[0]
    is_abs_path = root[0] == '/'
    clean = [str(s).strip('/') for s in args]
    if is_abs_path:
        clean[0] = '/' + clean[0]
    return '/'.join(clean)

def get_base_dir(saved_data, *args):
    network_root = get_dir_root(saved_data)
    return concat_path(network_root, '/mat/Fahimehb/bmtk_run', *args)


def get_output_dir(input_type, stim_type, model_type, saved_data, network_id, *args):
    """ Get dir containing runs for given params """
    base_dir = get_base_dir(saved_data)
    return concat_path(base_dir, '/outputs/', input_type, stim_type, model_type, network_id, *args)

def get_output_filename(amp, trial, freq=None):
    if not freq:
        return concat_path('amp' + format_amp(amp) + '_tr' + str(trial))
    else:
        return concat_path('amp' + format_amp(amp) + '_freq' + format_freq(freq) + '_tr' + str(trial))

def get_config_resolved_path(stim_type, out_folder, amp, freq=None):
        if stim_type == "dc":
            key = resolve_dc_key(amp)
        if stim_type == "sin":
            key = resolve_sin_key(amp, freq)

        return concat_path(out_folder, 'config_' + key + '.json')

def resolve_nodes_filename(network_id):
    parts = [network_id + '_nodes.h5']
    return  ''.join(parts)


def get_dir_name(network_id, amp, freq=None, trial=0):
    if freq is not None:
        return get_sin_dir_name(network_id, amp, freq, trial)
    else:
        return get_dc_dir_name(network_id, amp, trial)

### DC ###

def resolve_dc_key(amp):
    """ file/folder name code """
    parts = ['amp' + format_amp(amp)]
    return '_'.join(parts)

def get_dc_dir_name(network_id, amp, trial):
    return '_'.join([network_id , resolve_dc_key(amp), 'tr' + str(trial)])

### SIN ###

def resolve_sin_key(amp, freq):
    """ file/folder name code """
    parts = ['amp' + format_amp(amp), 'freq' + format_freq(freq)]
    return '_'.join(parts)

def get_sin_dir_name(network_id , amp, freq, trial):
    return '_'.join([network_id , resolve_sin_key(amp, freq), 'tr' + str(trial)])

### Table ###

def get_table_filename(network_id, amp, trial, freq=None):
    if freq is not None:
        return 'table_{}_amp{}_freq{}_tr{}.h5'.format(network_id, format_amp(amp), format_freq(freq), trial)
    else:
        return 'table_{}_amp{}_tr{}.h5'.format(network_id, format_amp(amp), trial)

def get_table_dir(input_type, stim_type, model_type, saved_data, *args):
        root_dir = get_dir_root(saved_data)
        return concat_path(root_dir, '/mat/Fahimehb/bmtk_run/result_tables/', input_type,
                              stim_type, model_type, *args)


#################################################
#
#     IO / data
#
#################################################


def get_electrode_xyz(electrode_pos_path):  # Ideally you would use the method bionet uses
    # mesh files are unnecessary for this study
    electrode_pos_df = pd.read_csv(electrode_pos_path, sep=' ')
    return electrode_pos_df.as_matrix(columns=['pos_x', 'pos_y', 'pos_z'])

def get_node_id(nodeh5, network_id):
    return nodeh5['nodes/'+network_id+'/node_id'].value

def get_node_type_id(nodeh5, network_id):
    return nodeh5['nodes/'+network_id+'/node_type_id'].value

def get_node_group_id(nodeh5, network_id):
    return nodeh5['nodes/'+network_id+'/node_group_id'].value

def get_node_group_index(nodeh5, network_id):
    return nodeh5['nodes/'+network_id+'/node_group_index'].value

def get_cell_xyz(nodeh5, network_id, cell_group_id):
    # print 'nodes/'+network_id+'/'+str(cell_group_id)+'/positions'
    return  nodeh5['nodes/'+network_id+'/'+str(cell_group_id)+'/positions'].value

def get_cell_spikes(spikesh5, gid):
    all_spike_gids = spikesh5['spikes/gids'].value
    all_spike_times = spikesh5['spikes/timestamps'].value
    df = pd.DataFrame(columns=["gid_spike", "gid_spike_time"])
    df['gid_spike'] = all_spike_gids
    df['gid_spike_time'] = all_spike_times
    return df[df['gid_spike'] == gid]['gid_spike_time'].tolist()

def get_vm_data(cvh5):
    gids = cvh5['mapping/gids'].value
    vm= cvh5['v/data'].value
    df = pd.DataFrame(vm, columns=gids)
    return df



def get_vext_data(cvh5):
    gids = cvh5['mapping/gids'].value
    vext= cvh5['vext/data'].value
    df = pd.DataFrame(vext, columns=gids)
    return df

def get_dt(conf):
    return conf['run']['dt']

def get_tstop(conf):
    return conf['run']['tstop']

def get_spikes_file(output):
    return h5.File(concat_path(output, 'spikes.h5'), 'r')


def get_cv_files(output):
    return h5.File(concat_path(output, 'cell_vars.h5'), 'r')

def get_nodes_files(network_id, run_dir):
    network_filename = resolve_nodes_filename(network_id)
    return h5.File(concat_path(run_dir, 'new_network' ,network_filename), 'r')

def get_json_from_file(path):
    with open(path, 'r') as f:
        return json.load(f)



#################################################
#
#     Run aggregate data output file
#
#################################################

dc_cols = ['trial', 'gid', 'node_type_id', 'x', 'y', 'z', 'distance', 'amp', 'spikes']
sin_cols=['fq']
spike_phase_analysis_cols=['spike_threshold_t', 'spike_phase']


def resolve_additional_cols(include_sin, include_spike_phase):
    """ Find all the additional columns for the table"""
    cols = []
    if include_sin:
        cols= cols + sin_cols
    # if include_spike_phase:
    #     cols = cols + spike_analysis_cols
    return cols


def build_df(additional_cols=[]):
    """ Wrapped df creation to give place to explicitly declare column types """
    df = pd.DataFrame(columns=(dc_cols + additional_cols))
    # pd doesn't do a great job of identifying ints
    df['trial'] = df['trial'].astype(int)
    return df

def resolve_run_id(gid, amp, freq=None):

    if freq is not None:
        stringified = map(str, [gid, format_amp(amp), freq])
    else:
        stringified = map(str, [gid, format_amp(amp)])
    return '_'.join(stringified)

def resolve_output_aggregates(network_id, include_sin,  amp, trial):
    """ Find all the outputs for the table"""
    if  not include_sin:
        return get_dir_name(network_id= network_id, amp=amp, trial=trial)
    if include_sin:
        return get_dir_name(network_id = network_id, amp=amp, freq=999, trial=trial) .replace('999', '*')



def write_table_h5(fpath, df, attrs=None):
   """ dataframe to h5 """
   # df.sort_values('gid', inplace=True)
   with h5.File(fpath, 'w') as f5:
       f5.create_dataset('ids', data=map(str, df.index))
       for col in df.columns:
           if col not in ['spikes', 'spike_threshold_t', 'spike_phase']: # cuz it will explode if it is type 'object'
               f5.create_dataset(col, data=df[col].astype(float))

       spike_grp = f5.create_group('spikes')
       for (rid, spike_times) in df.spikes.iteritems():
           spike_grp.create_dataset(str(rid), maxshape=(None,), chunks=True, data=spike_times)


       if set(['spike_threshold_t', 'spike_phase']).issubset(df.columns):
           spike_threshold_t_grp = f5.create_group('spike_threshold_t')
           spike_phase_grp = f5.create_group('spike_phase')

           for (rid, s_t_t) in df.spike_threshold_t.iteritems():
               spike_threshold_t_grp.create_dataset(str(rid), maxshape=(None,), chunks=True, data=s_t_t)
       #
           for (rid, s_p) in df.spike_phase.iteritems():
               spike_phase_grp.create_dataset(str(rid), maxshape=(None,), chunks=True, data=s_p)

       if attrs is not None:
           for k,v in attrs.iteritems():
               f5.attrs[k] = v





def read_table_h5(fpath):
    """ h5 to dataframe """

    with h5.File(fpath, 'r') as f5:
        extra_cols = []
        if 'has_sin' in f5.attrs and f5.attrs['has_sin']:
            extra_cols = extra_cols + sin_cols

        if 'has_spike_phase_analysis' in f5.attrs and f5.attrs['has_spike_phase_analysis']:
             extra_cols = extra_cols + spike_phase_analysis_cols

        table = build_df(extra_cols)
        un_touched_data_cols = (dc_cols + extra_cols)

        ids = f5['ids'].values
        spike_data = f5['spikes']

        if set(['spike_threshold_t', 'spike_phase']).issubset(un_touched_data_cols):
             spike_threshold_t_data = f5['spike_threshold_t']
             spike_phase_data = f5['spike_phase']

        for i, rid in enumerate(ids):  # i is key for f5 dsets, rid is key for spikes & df
            data_cols = (dc_cols + extra_cols)
            spike_index = data_cols.index('spikes')

            if set(['spike_threshold_t', 'spike_phase']).issubset(data_cols):
                spike_threshold_t_index = data_cols.index('spike_threshold_t')
                spike_phase_index = data_cols.index('spike_phase')

            data_cols.pop(spike_index)
            if set(['spike_threshold_t', 'spike_phase']).issubset(data_cols):
                new_spike_threshold_t_index = data_cols.index('spike_threshold_t')
                data_cols.pop(new_spike_threshold_t_index)
                new_spike_phase_index = data_cols.index('spike_phase')
                data_cols.pop(new_spike_phase_index)

            data = [f5[c][i] for c in data_cols]
            # place spikes in correct position
            data.insert(spike_index, spike_data[rid].value)
            if set(['spike_threshold_t', 'spike_phase']).issubset(un_touched_data_cols):
                data.insert(spike_threshold_t_index, spike_threshold_t_data[rid].value)
                data.insert(spike_phase_index, spike_phase_data[rid].value)

            table.loc[rid] = data

    return table


def read_network_tables(network_id, amp_range, input_type, stim_type, model_type, trial, freq , saved_data=False, data_dir=None):
    """ Read h5 files for a set of amplitudes """

    print "Fetching data..."
    data_dir = get_table_dir(input_type, stim_type, model_type, saved_data=saved_data) if data_dir is None else data_dir
    paths = [concat_path(data_dir, get_table_filename(network_id, amp, trial, freq=freq)) for amp in amp_range]
    t = build_df()  # do this for code analysis
    t = t.append([read_table_h5(p) for p in paths])
    # t['num_spikes'] = t.apply(lambda row: len(row['spikes']), axis=1)
    # t['num_true_spikes'] = np.where(t['num_spikes'] == 1, 0, t['num_spikes'])

    # print "Done"
    return t

#################################################
#
#     ????
#
#################################################
def get_spike_h5_file(output):
    return concat_path(output, "spikes.h5")

def get_cellvars_h5_file(output):
    return concat_path(output, "cell_vars.h5")

def get_spike_time(spikes_h5_file):
    return  spikes_h5_file['spikes/timestamps'].value

def get_gids(spikes_h5_file):
    return spikes_h5_file['spikes/gids'].value

