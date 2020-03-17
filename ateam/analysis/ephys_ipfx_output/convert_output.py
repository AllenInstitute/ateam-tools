import json
import numpy as np
import pandas as pd
import ipfx.bin.lims_queries as lq
import allensdk.core.json_utilities as ju

cell_features = [
    'vrest',
    'ri',
    'sag',
    'tau',
    'f_i_curve_slope',
    # 'adaptation',
    # 'latency',
    # 'avg_isi',
    # 'vm_for_sag',
    # 'has_burst',
    # 'has_pause',
    # 'has_delay',
]
# from cell_features.long_squares
ls_features =[
    "input_resistance",
    "tau",
    "v_baseline",
    "sag",
    "rheobase_i",
    "fi_fit_slope",
]
# # all available:
# sweep_feature_list = [
#     "first_isi",
#     "avg_rate",
#     "isi_cv",
#     "latency",
#     "median_isi",
#     "adapt",
# ]
hero_sweep_features = [
    # 'adapt',
    'avg_rate',
    'first_isi',
    # 'mean_isi',
    # 'latency',
]
rheo_sweep_features = [
    'latency',
    # 'first_isi',
]
mean_sweep_features = [
    'adapt',
]
spike_features = [
    'upstroke_downstroke_ratio',
    'threshold_v',
    'peak_v',
    'fast_trough_v',
    'trough_v',
    # 'threshold_i',
    # include all troughs?
    # 'slow_trough_v',
]
ls_spike_features = [
    'upstroke_downstroke_ratio',
    'threshold_v',
    'peak_v',
    'fast_trough_v',
    'trough_v',
    # not in cell record
    'width',
    'upstroke',
    'downstroke',
]
invert = ["first_isi"]
spike_threshold_shift = ["trough_v", "fast_trough_v", "peak_v"]

sweep_types = [
    '_long_square',
    # '_short_square', 
    # '_ramp', 
    ]

def extract_local_pipeline_output(output_json):
    output = ju.read(output_json)
    record = {}
    cell_state = output.get('qc', {}).get('cell_state')
    if cell_state is not None:
        record['failed_qc'] = cell_state.get('failed_qc', False)
        record['fail_tags'] = '; '.join(cell_state.get('fail_tags'))

    fx_dict = output.get('feature_extraction')
    if fx_dict is not None:
        record.update(extract_fx_output(fx_dict))
    return record

def extract_fx_output(fx_dict, v2=False):
    record = {}
    cell_state = fx_dict.get('cell_state')
    if cell_state is not None:
        record['failed_fx'] = cell_state.get('failed_fx', False)
        fail_message = cell_state.get('fail_fx_message')
        # TODO: can remove list option after offpipeline rerun
        record['fail_fx_message'] = '; '.join(fail_message) if isinstance(fail_message, list) else fail_message

    # cell_record = fx_dict.get('cell_record')
    # if cell_record is not None:
    #     record.update({key: cell_record.get(key) for key in cell_features})

    # already averaged first spike features
    # for feature in spike_features:
    #     for sweep_type in sweep_types:
    #         key = feature+sweep_type
    #         record[key] = cell_record.get(key)
    if v2:
        cell_features = fx_dict["specimens"][0].get('cell_ephys_features', {})
    else:
        cell_features = fx_dict.get('cell_features', {})
    # if cell_features["ramps"]:
    #     ramp_mean_spike_0 = cell_features["ramps"]["mean_spike_0"]
    #     add_features_to_record(spike_features, ramp_mean_spike_0, record, suffix="_ramp")

    # if cell_features["short_squares"]:
    #     sq_mean_spike_0 = cell_features["short_squares"]["mean_spike_0"]
    #     add_features_to_record(spike_features, sq_mean_spike_0, record, suffix="_short_square")

    long_squares = cell_features.get('long_squares')
    if long_squares is not None:
        add_features_to_record(ls_features, long_squares, record)

        sweep = long_squares.get('rheobase_sweep',{})
        add_features_to_record(rheo_sweep_features, sweep, record, suffix='_rheo')
        ls_rheo_spike_0 = sweep["spikes"][0]
        add_features_to_record(ls_spike_features, ls_rheo_spike_0, record, suffix="_long_square")

        sweep = long_squares.get('hero_sweep',{})
        add_features_to_record(hero_sweep_features, sweep, record, suffix='_hero')


        sweeps = long_squares.get('spiking_sweeps',{})
        if sweeps is not None:
            for feature in mean_sweep_features:
                key = feature+'_mean'
                feat_list = [sweep[feature] for sweep in sweeps if feature in sweep]
                record[key] = np.mean([x for x in feat_list if x is not None])

    offset_feature_values(spike_threshold_shift, record, "threshold_v")
    invert_feature_values(invert, record)
    return record

def offset_feature_values(features, record, relative_to):
    for feature in features:
        matches = [x for x in record if x.startswith(feature)]
        for match in matches:
            suffix = match[len(feature):]
            val = record.pop(match)
            record[match+"_rel"] = (val - record[relative_to+suffix]) if val is not None else None

def invert_feature_values(features, record):
    for feature in features:
        matches = [x for x in record if x.startswith(feature)]
        for match in matches:
            suffix = match[len(feature):]
            val = record.pop(match)
            record[match+"_inv"] = 1/val if val is not None else None

def add_features_to_record(features, feature_data, record, suffix=""):
    record.update({feature+suffix: feature_data.get(feature) for feature in features})

def compile_lims_results(specimen_ids):
    records = []
    for cell in specimen_ids:
        path = lq.get_fx_output_json(cell)
        if path.startswith('/'):
            record = extract_fx_output(ju.read(path), v2=("V2" in path))
            record["specimen_id"] = cell
            records.append(record)
    ephys_df = pd.DataFrame.from_records(records, index="specimen_id")
    return ephys_df
