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
    'avg_isi',
    'f_i_curve_slope',
    'adaptation',
    'latency',
    # 'vm_for_sag',
    # 'has_burst',
    # 'has_pause',
    # 'has_delay',
]

sweep_features = [
    'adapt'
]

spike_features = [
    'upstroke_downstroke_ratio',
    'peak_v',
    'trough_v',
    'fast_trough_v',
    'slow_trough_v',
    'threshold_v',
    'threshold_i',
]

sweep_types = ['_ramp', '_short_square', '_long_square']

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

def extract_fx_output(fx_dict):
    record = {}
    record['failed_fx'] = fx_dict['cell_state'].get('failed_fx', False)
    fail_message = fx_dict['cell_state'].get('fail_fx_message')
    record['fail_fx_message'] = '; '.join(fail_message) if isinstance(fail_message, list) else fail_message

    cell_record = fx_dict.get('cell_record')
    if cell_record is not None:
        record.update({key: cell_record.get(key) for key in cell_features})

        # already averaged first spike features
        for feature in spike_features:
            for sweep_type in sweep_types:
                key = feature+sweep_type
                record[key] = cell_record.get(key)

        cf = fx_dict.get('cell_features', {})
        long_squares = cf.get('long_squares')
        if long_squares is not None:
            sweeps = long_squares.get('spiking_sweeps')
            if sweeps is not None:
                for feature in sweep_features:
                    key = feature+'_mean'
                    feat_list = [sweep[feature] for sweep in sweeps if feature in sweep]
                    record[key] = np.mean([x for x in feat_list if x is not None])

    return record
    
def compile_lims_results(specimen_ids):
    records = []
    for cell in specimen_ids:
        path = lq.get_fx_output_json(cell)
        if path.startswith('/'):
            record = extract_fx_output(ju.read(path))
            record["specimen_id"] = cell
            records.append(record)
    ephys_df = pd.DataFrame.from_records(records, index="specimen_id")
    return ephys_df
