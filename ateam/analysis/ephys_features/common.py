import numpy as np
import pandas as pd
from itertools import chain
import warnings
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor

# suppress output from AllenSDK deprecated use of dataset.value
from h5py.h5py_warnings import H5pyDeprecationWarning
warnings.filterwarnings("ignore", category=H5pyDeprecationWarning)

def spike_times(v, i, t, start=None, end=None):
    sweepex = EphysSweepFeatureExtractor(t=t, v=v, i=i, start=start, end=end)
    try:
        sweepex.process_spikes()
        return sweepex.spike_feature("peak_t")
    except ValueError:
        warnings.warn('Error processing spikes')
        return []

def stim_props(v, i, t):
    didt = np.diff(i)/(t[1]-t[0])
    ichange = np.flatnonzero(np.abs(didt) > 1e-3) # 1 pA/ms
    # past initial test pulse
    ichange = ichange[t[ichange] > 0.1]
    # nonconsecutive
    # ichange = ichange[np.flatnonzero(np.diff(ichange) > 1) + 1]
    tstim = t[ichange]
    duration = tstim[1] - tstim[0]
    amp = np.mean(i[ichange[0]+1 : ichange[1]-1])
    if (len(tstim) > 2) or (abs(duration - 1) > 0.01):
        warnings.warn('Multiple stimulus pulses')
    return tstim, amp

def grouped_agg_flat(groupby, agg_functions=['mean']):
    df = groupby.agg(agg_functions)
    df.columns = ['_'.join(col).rstrip('_') for col in df.columns.values]
    return df

def perithreshold(sweep_results_df):
    current = "current"
    df = sweep_results_df.sort_values(current, ascending=False).groupby("cell_id").nth(0)
    df.columns = [col+'_perithresh' for col in df.columns.values]
    return df

def process_sweep_results(sweep_results_df, save_path=None):
    agg = grouped_agg_flat(sweep_results_df.groupby("cell_id"), agg_functions=['mean','max','min','median'])
    perithresh = perithreshold(sweep_results_df)
    cell_results_df = pd.concat([agg, perithresh], axis=1)
    if save_path:
        sweep_results_df.to_csv(save_path + "_sweeps.csv")
        cell_results_df.to_csv(save_path + "_cells.csv")
    return cell_results_df

def run_ephys_analysis(cell_function, specimen_ids, processes=1, save_path=None):
    try:
        from multiprocessing import Pool
        p = Pool(processes=processes)
        results = p.map(cell_function, specimen_ids)
    except ImportError:
        results = map(cell_function, specimen_ids)

    sweep_results_df = pd.DataFrame.from_records(chain(*results))
    ca_df = process_sweep_results(sweep_results_df, save_path=save_path)
    return ca_df