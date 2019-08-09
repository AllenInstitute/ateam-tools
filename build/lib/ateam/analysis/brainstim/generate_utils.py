######################################################
# Authors: Taylor Connington
# Date created: 9/1/2017
######################################################

import shutil
import copy
from reobase_utils import *
import itertools

"""
Set of functions strictly for generating batches of config files. The separation from reobase_utils is a bit arbitrary,
but the functions are long as stand on their own. It should be separate from the script however, because then this 
logic is reusable by another script if required (e.g. a parallel batch run script or something)
"""

def generate_config_set(config_base, bconf, confs_folder='confs', verbose=False):
    """ Parse batch config file and enumerate/use all run combinations; run validation checks """

    validate_bconf(bconf) # try to avoid wonky stuff

    cell = bconf['cell_gid']
    el_range = bconf['el_range']
    els = range(el_range[0], el_range[1])
    considerations = [els, bconf['amps']]
    stim_type = bconf['stim_type']
    model_type = bconf['model_type']
    input_type = bconf['input_type']
    extra_waveform_type = bconf['stim_type'].split('_')[0]

    # key function used for consistent naming
    key_fn = {
        'dc': resolve_dc_key,
        'sin': resolve_sin_key,
        'sin_dc': resolve_sin_dc_key,
    }.get(stim_type)
    # dir function used for consistent naming
    dir_fn = {
        'dc': get_dc_dir_name,
        'sin': get_sin_dir_name,
        'sin_dc': get_sin_dc_dir_name,
    }.get(stim_type)

    # only add freqs if it's there
    if 'freqs' in bconf:
        considerations += [bconf['freqs']]
    if 'ic_amps' in bconf:
        considerations += [bconf['ic_amps']]

    # empty lists will cause empty product below
    for arr in considerations:
        if len(arr) == 0:
            raise ValueError('Cannot pass empty array as value set in for batch config')

    for combo in itertools.product(*considerations):
        el = combo[0]
        key = key_fn(*combo)
        run_folder = dir_fn( *(combo + (bconf["trial"],)) )
        extra_waveform_props = {
            "shape": extra_waveform_type,
            "amp": combo[1],
        }

        if extra_waveform_type == 'sin':
            extra_waveform_props["freq"] = combo[2]

        if input_type == 'extrastim_intrastim':
            Iclamp_props = {
                "amp": combo[3],
            }

        conf_name = 'config_{}.json'.format(key)


        conf_data = generate_config(config_base, conf_name, confs_folder,
                                    el, cell, run_folder, extra_waveform_props, stim_type, model_type, input_type)
        if input_type == 'extrastim_intrastim':
            conf_data = generate_config(config_base, conf_name, confs_folder,
                                        el, cell, run_folder, extra_waveform_props, stim_type, model_type, input_type, Iclamp_props)

        if verbose:
            print 'KEY: ', key, '   -----   OUTPUT DIR: ', conf_data["manifest"]["$OUTPUT_DIR"]


def set_config(conf_data, el, cell_gid, run_folder, extra_waveform_props, stim_type, model_type, input_type, Iclamp_props=None):
    """ Lots of config settings are cell specific or input specific--all changes to the conf file should go here """
    cell_gid = str(cell_gid)
    # run_folder = get_dc_dir_name(el, amp, trial)

    conf_data["extracellular_stimelectrode"]["position"] = "$STIM_DIR/" + cell_gid + "_" + format_el(el) + ".csv"
    conf_data["extracellular_stimelectrode"]["waveform"].update(extra_waveform_props)
    if Iclamp_props is not None:
        conf_data["iclamp"].update(Iclamp_props)
    # Note: output dir doesn't include current sign
    outdir_root = conf_data["manifest"]["$OUTPUT_DIR"]
    conf_data["manifest"]["$OUTPUT_DIR"] = concat_path(outdir_root, input_type, stim_type, model_type, cell_gid, run_folder)
    conf_data["manifest"]["$STIM_DIR"]   = concat_path("$BASE_DIR/stimulation", cell_gid)
    # single cell definition and per-cell connection functions
    conf_data["internal"]["cells"]     = "$NETWORK_DIR/{}_cell.csv".format(cell_gid)
    conf_data["internal"]["con_types"] = "$NETWORK_DIR/func.csv".format(cell_gid)

    return conf_data

def generate_config(base, filename, out_dir, *args, **kwargs):
    """ Per file logic -- modify base config data and then store in new file """
    with open(base, 'r') as fp:
        base_data = json.load(fp)

    data = set_config(copy.deepcopy(base_data), *args, **kwargs)

    with open(out_dir + '/' + filename, 'w') as fp:
        json.dump(data, fp, indent=4, separators=(',', ': ')) # print pretty

    return data

def validate_bconf(bconf):

    must_haves = {'el_range':list, 'cell_gid':None, 'amps':list, 'trial':int, 'stim_type':str, 'model_type':str, 'input_type': str}

    for name, t in must_haves.iteritems():
        if name not in bconf:
            raise ValueError("You are required to include '{}' in the batch configuration".format(name))

        if name == 'cell_gid':
            if not bconf[name].isdigit():
                raise ValueError("Invalid cell_gid {}".format(bconf[name]))
        elif type(bconf[name]) != t:
            raise ValueError("Invalid type for '{}' ({})".format(name, type(bconf[name])))

    # input_type = bconf['input_type']
    extra_waveform_type = bconf['stim_type'].split('_')[0]
    input_type = bconf['input_type']

    if bconf['input_type'] != 'extrastim':
        intra_wavefrom_type = bconf['stim_type'].split('_')[1]
        print "Intrastim waveform is: ", intra_wavefrom_type

    print "Extrastim waveform is: ", extra_waveform_type

    if extra_waveform_type == 'dc':
        if 'freqs' in bconf:
            raise ValueError("Cannot include 'freqs' key in batch config when using DC input")
    elif extra_waveform_type == 'sin':
        if 'freqs' not in bconf:
            raise ValueError("Must include 'freqs' list in batch config when using sin input")
    else:
        print 'Cannot validate unknown extra_waveform type {} --I suggest implementing some common sense checks'.format(extra_waveform_type)


    if input_type == 'extrastim':
        if 'ic_amps' in bconf:
            raise ValueError("Cannot include 'ic_amps' key in batch config when using only extrastim input")
    elif input_type != 'extrastim':
        if 'ic_amps' not in bconf:
            raise ValueError("Must include 'ic_amps' list in batch config when using intrastim input")


# "Reset" a folder - would make a good general util
def prep_confs_folder(confs_folder):
    if os.path.isdir(confs_folder):
        # clear it
        shutil.rmtree(confs_folder)

    os.mkdir(confs_folder) # placed in current dir

