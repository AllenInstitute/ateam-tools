import pandas as pd 
import numpy as np
from ateam.data.convert_ipfx_output_full import add_features_to_record, get_cell_features_lims

mean_sweep_features = ['adapt',  'isi_cv', 'threshold_v', 'peak_v', 'trough_v', 'upstroke', 'downstroke', 'width']
interp_features = [ 'latency', 'avg_rate',]
subthreshold_features = ['v_baseline', 'sag', 'input_resistance']

def get_interpolated_long_square_features(cell_features, amplitudes):
    record = {}
    long_squares_analysis = cell_features.get("long_squares")
    if long_squares_analysis is None:
        return record
    add_features_to_record(subthreshold_features, long_squares_analysis, record)

    sweeps = long_squares_analysis.get('spiking_sweeps')
    for sweep in sweeps:
        sweep.update(sweep.pop("spikes")[0])
    sweeps_df = pd.DataFrame.from_records(sweeps)
    record.update(interpolate_feature_value(sweeps_df, interp_features, amplitudes))
    for feature in mean_sweep_features:
        record[feature] = sweeps_df[feature].mean()
    
    return record

def interpolate_feature_value(sweeps_df, features, amplitudes):
    amp_var = "stim_amp"
    sweeps_df.sort_values(amp_var, inplace=True)
    n_amp = len(amplitudes)
    record = {}
    for feature in features:
        values = np.interp(amplitudes, sweeps_df[amp_var], sweeps_df[feature], left=np.nan, right=np.nan)
        record.update({f"{feature}_{i}": values[i] for i in range(n_amp)})
    return record

def compile_lims_results(specimen_ids):
    records = []
    for cell in specimen_ids:
        cell_features = get_cell_features_lims(cell)
        if cell_features:
            record = get_interpolated_long_square_features(cell_features)
            record["specimen_id"] = cell
            records.append(record)
    ephys_df = pd.DataFrame.from_records(records, index="specimen_id")
    return ephys_df