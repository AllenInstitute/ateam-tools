import numpy as np
import pandas as pd
import argschema as ags
import glob
import os.path
import ipfx.lims_queries as lq
# from ipfx.stimulus_protocol_analysis import LongSquareAnalysis
import allensdk.core.json_utilities as ju
import logging

class CompileOutputParameters(ags.ArgSchema):
    ids = ags.fields.List( ags.fields.Integer,
        description="List of specimen IDs to process",
        cli_as_single_argument=True,
        default=None, allow_none=True
    )
    input_file = ags.fields.InputFile(
        description=("Input file of specimen IDs (one per line)"
                     "- optional if LIMS is source"),
        default=None, allow_none=True
    )
    project = ags.fields.String(
        description="Project code used for LIMS query, e.g. hIVSCC-MET",
        default="hIVSCC-MET",
        allow_none=True
    )
    cell_count_limit = ags.fields.Integer(
        description="Limit to number of cells evaluated",
        default=1000
    )
    output_file = ags.fields.OutputFile(
        description="output file path",
        default="compiled_output_ipfx_lims.csv",
        allow_none=True
    )

# from cell_features.long_squares
ls_features = [
    "input_resistance",
    "tau",
    "v_baseline",
    "sag",
    "rheobase_i",
    "fi_fit_slope",
]
hero_sweep_features = [
    'adapt',
    'avg_rate',
    'latency',
    'first_isi',
    'mean_isi',
    "median_isi",
    "isi_cv",
]
rheo_sweep_features = [
    'latency',
    'avg_rate',
]
mean_sweep_features = [
    'adapt',
]
ss_spike_features = [
    'upstroke_downstroke_ratio',
    'threshold_v',
    'peak_v',
    'fast_trough_v',
]
ramp_spike_features = [
    'upstroke_downstroke_ratio',
    'threshold_v',
    'peak_v',
    'fast_trough_v',
    'trough_v',
    'threshold_i',
]
ls_spike_features = [
    'upstroke_downstroke_ratio',
    'threshold_v',
    'peak_v',
    # include all troughs?
    'fast_trough_v',
    'trough_v',
    # not in cell record
    'width',
    'upstroke',
    'downstroke',
]
spike_adapt_features = [
    'isi',
    'width',
    'upstroke',
    'downstroke',
    'threshold_v',
    'fast_trough_v',
]
invert = ["first_isi"]
spike_threshold_shift = ["trough_v", "fast_trough_v", "peak_v"]


def extract_local_pipeline_output(output_json, extra=[]):
    output = ju.read(output_json)
    record = {}
    cell_state = output.get('qc', {}).get('cell_state')
    if cell_state is not None:
        record['failed_qc'] = cell_state.get('failed_qc', False)
        record['fail_tags'] = '; '.join(cell_state.get('fail_tags'))

    fx_dict = output.get('feature_extraction')
    if fx_dict is not None:
        record.update(extract_fx_output(fx_dict, extra=extra))
    return record

def get_cell_features_lims(specimen_id):
    path = get_fx_output_json(specimen_id)
    if not path.startswith('/'):
        return None
    v2 = ("V2" in path) or ("DATAFIX" in path)
    fx_dict = ju.read(path)
    if v2:
        cell_features = fx_dict["specimens"][0].get('cell_ephys_features', {})
    else:
        cell_features = fx_dict.get('cell_features', {})
    return cell_features

def extract_fx_output(fx_dict, v2=False, extra=['60pa','5spikes']):
    record = {}
    cell_state = fx_dict.get('cell_state')
    if cell_state is not None:
        # record['failed_fx'] = cell_state.get('failed_fx', False)
        record['fail_fx_message'] = cell_state.get('fail_fx_message')

    if v2:
        cell_features = fx_dict["specimens"][0].get('cell_ephys_features', {})
    else:
        cell_features = fx_dict.get('cell_features', {})

    ramps = cell_features.get('ramps')
    if ramps is not None:
        mean_spike_0 = ramps["mean_spike_0"]
        add_features_to_record(ramp_spike_features, mean_spike_0, record, suffix="_ramp")

    short_squares = cell_features.get('short_squares')
    if short_squares is not None:
        mean_spike_0 = short_squares["mean_spike_0"]
        add_features_to_record(ss_spike_features, mean_spike_0, record, suffix="_short_square")

    # spiking_sweep_features_df = pd.DataFrame.from_records(sweeps)
    # if '60pa' in extra:
    #     sweep = LongSquareAnalysis.find_hero_sweep(long_squares_analysis['rheobase_i'], spiking_sweep_features_df, min_offset=59)
    #     if sweep is not None:
    #         add_features_to_record(hero_sweep_features, sweep, record, suffix='_60pa')
    #         add_features_to_record(ls_spike_features, sweep["spikes"][0], record, suffix="_60pa")
    #         add_feature_ratios_to_record(spike_adapt_features, sweep["spikes"], record, suffix="_60pa")

    # if '5spikes' in extra:
    #     spike_sets = [sweep["spikes"] for sweep in sweeps]
    #     sweep = find_spiking_sweep_by_min_spikes(spiking_sweep_features_df, spike_sets)
    #     if sweep is not None:
    #         add_features_to_record(hero_sweep_features, sweep, record, suffix='_5spikes')
    #         add_features_to_record(ls_spike_features, sweep["spikes"][0], record, suffix="_5spikes")
    #         add_feature_ratios_to_record(spike_adapt_features, sweep["spikes"], record, suffix="_5spikes")
            
    offset_feature_values(spike_threshold_shift, record, "threshold_v")
    invert_feature_values(invert, record)

    long_squares_analysis = cell_features.get('long_squares')
    if long_squares_analysis is not None:
        record.update(get_complete_long_square_features(long_squares_analysis))
    
    return record

def get_complete_long_square_features(long_squares_analysis):
    record = {}
    add_features_to_record(ls_features, long_squares_analysis, record)

    sweep = long_squares_analysis.get('rheobase_sweep',{})
    add_features_to_record(rheo_sweep_features, sweep, record, suffix='_rheo')
    add_features_to_record(ls_spike_features, sweep["spikes"][0], record, suffix="_rheo")

    sweep = long_squares_analysis.get('hero_sweep',{})
    add_features_to_record(hero_sweep_features, sweep, record, suffix='_hero')
    add_features_to_record(ls_spike_features, sweep["spikes"][0], record, suffix="_hero")
    add_feature_ratios_to_record(spike_adapt_features, sweep["spikes"], record, suffix="_hero")

    sweeps = long_squares_analysis.get('spiking_sweeps',{})
    # TODO: work on dataframe / reuse code
    for feature in mean_sweep_features:
        key = feature+'_mean'
        feat_list = [sweep[feature] for sweep in sweeps if feature in sweep]
        record[key] = np.mean([x for x in feat_list if x is not None])
    
    
    offset_feature_values(spike_threshold_shift, record, "threshold_v")
    invert_feature_values(invert, record)


def offset_feature_values(features, record, relative_to):
    for feature in features:
        matches = [x for x in record if x.startswith(feature)]
        for match in matches:
            suffix = match[len(feature):]
            val = record.pop(match)
            feature_short = feature[:-2] #drop the "_v"
            record[feature_short + "_deltav" + suffix] = (val - record[relative_to+suffix]) if val is not None else None

def invert_feature_values(features, record):
    for feature in features:
        matches = [x for x in record if x.startswith(feature)]
        for match in matches:
            suffix = match[len(feature):]
            val = record.pop(match)
            record[feature + "_inv" + suffix] = 1/val if val is not None else None

def add_features_to_record(features, feature_data, record, suffix=""):
    record.update({feature+suffix: feature_data.get(feature) for feature in features})

def add_feature_ratios_to_record(features, dict_list, record, suffix=""):
    suffix = "_adapt_ratio" + suffix
    nspikes = len(dict_list)
    for i in range(nspikes-1):
        dict_list[i]['isi'] = dict_list[i+1]['peak_t'] - dict_list[i]['peak_t']
    if nspikes<4:
        return
    for feature in features:
        last = dict_list[-2].get(feature)
        # if last is None and nspikes >= 3:
        #     last = dict_list[-2].get(feature)
        if last is not None:
            record.update({feature+suffix: last/dict_list[0].get(feature)})

def find_spiking_sweep_by_min_spikes(spiking_features, spikes_set, min_spikes=5):
    num_spikes = np.array([len(spikes) for spikes in spikes_set])
    # spiking_features['spikes'] = spikes_set
    spiking_features = spiking_features.loc[num_spikes >= min_spikes].sort_values("stim_amp")
    spiking_features_depolarized = spiking_features[spiking_features["stim_amp"] > 0]

    if spiking_features_depolarized.empty:
        logging.info("Cannot find sweep with >{min_spikes} spikes.")
        return None
    else:
        return spiking_features_depolarized.iloc[0]

def compile_lims_results(specimen_ids, extra=['60pa','5spikes']):
    records = []
    for cell in specimen_ids:
        path = get_fx_output_json(cell)
        if path.startswith('/'):
            record = extract_fx_output(ju.read(path), v2=("V2" in path) or ("DATAFIX" in path), extra=extra)
            record["specimen_id"] = cell
            records.append(record)
    ephys_df = pd.DataFrame.from_records(records, index="specimen_id")
    return ephys_df

def get_specimen_ids(ids=None, input_file=None, project="T301", include_failed_cells=False, cell_count_limit=float('inf')):
    if ids is not None:
        specimen_ids = ids
    elif input_file is not None: 
        with open(module.args["input_file"], "r") as f:
            ids = [int(line.strip("\n")) for line in f]
        module.args.pop('ids')
    else:
        specimen_ids = lq.project_specimen_ids(
            project, passed_only=not include_failed_cells)
    if len(specimen_ids) > cell_count_limit:
        specimen_ids = specimen_ids[:cell_count_limit]
    logging.info(
        "Number of specimens to process: {:d}".format(len(specimen_ids)))
    return specimen_ids

def get_fx_output_json(specimen_id):
    """
    Find in LIMS the full path to the json output of the feature extraction module
    If more than one file exists, then chose the latest version

    Parameters
    ----------
    specimen_id

    Returns
    -------
    file_path: string
    """
    NO_SPECIMEN = "No_specimen_in_LIMS"
    NO_OUTPUT_FILE = "No_feature_extraction_output"
    
    sql = """
    select err.storage_directory, err.id
    from specimens sp
    join ephys_roi_results err on err.id = sp.ephys_roi_result_id
    where sp.id = %d
    """ % specimen_id

    res = lq.query(sql)
    if res:
        err_dir = res[0]["storage_directory"]

        file_list = glob.glob(os.path.join(err_dir, '*EPHYS_FEATURE_EXTRACTION_*_output.json'))
        if file_list:
            latest_file = max(file_list, key=os.path.getctime)   # get the most recent file
            return latest_file
        else:
            return NO_OUTPUT_FILE
    else:
        return NO_SPECIMEN

def main(ids=None, input_file=None, project="T301", include_failed_cells=False, cell_count_limit=float('inf'), output_file=None, **kwargs):
    specimen_ids = get_specimen_ids(ids, input_file, project, include_failed_cells, cell_count_limit)
    compile_lims_results(specimen_ids).to_csv(output_file)

if __name__ == "__main__":
    module = ags.ArgSchemaParser(schema_type=CompileOutputParameters)
    main(**module.args)