import numpy as np
from ipfx.sweep import Sweep, SweepSet
import ipfx.stim_features as stf
from ipfx import feature_extractor as fx
from ipfx.stimulus_protocol_analysis import StimulusProtocolAnalysis
from ateam.analysis.bmtk.cell_vars import get_cellvar_report
import pandas as pd

def get_voltage_features_config(config_file, start=0, end=None):
    var_report = get_cellvar_report(config_file)
    gids = var_report.gids

    records = []
    analysis = SingleStepAnalysis(start=start/1000, end=end/1000)
    for gid in gids:
        # , time_window=[start, end]
        v, t = var_report.data_timeseries(gid, 'v', compartments='origin')
        i = np.zeros_like(v)
        sweep_set = sweep_set_for_model(t/1000., v, i)
        analysis.analyze(sweep_set)
        features = analysis.sweep_features()
        features.update(analysis.spike_ratio_features(['width']))
        records.append(features)
    return pd.DataFrame.from_records(records, index=gids)

def sweep_set_for_model(t, v, i):
    """Generate a SweepSet object based on a single model sweep
    Parameters
    ----------
    t: array
        Time data (seconds)
    v: array
        Voltage data
    i: array
        Current stimulus data
    Returns
    -------
    SweepSet containing one Sweep
    """
    sampling_rate = 1 / (t[1] - t[0])
    sweep = Sweep(t=t,
                  v=v,
                  i=i,
                  sampling_rate=sampling_rate,
                  sweep_number=None,
                  clamp_mode="CurrentClamp",
                  epochs=None,
                  )
    return SweepSet([sweep])


class SingleStepAnalysis(StimulusProtocolAnalysis):
    """ Analysis of responses to step current stimuluation
    Parameters
    ----------
    start: float
        Start time of stimulus interval (seconds)
    end: float
        End time of stimulus interval (seconds)
    """
    def __init__(self, start=None, end=None):
        spx = fx.SpikeFeatureExtractor(start=start, end=end)
        sptx = fx.SpikeTrainFeatureExtractor(start, end,
            stim_amp_fn=stf._step_stim_amp)
        super(SingleStepAnalysis, self).__init__(spx, sptx)

    def analyze(self, sweep_set, exclude_clipped=False):
        """ Analyze spike and sweep features
        Parameters
        ----------
        sweep_set: SweepSet
        exclude_clipped: bool (optional, default=False)
            Whether to exclude clipped spikes from sweep-level features
        """
        extra_sweep_features = ["stim_amp", "v_baseline"]
        self.analyze_basic_features(sweep_set,
            extra_sweep_features=extra_sweep_features,
            exclude_clipped=exclude_clipped)

        # Analyze additional spike-level features
        # for sd in self._spikes_set:
        #     if sd.shape[0] >= 2:
        #         sd["slow_trough_delta_v"] = _slow_trough_delta_v(
        #             sd["fast_trough_v"].values, sd["slow_trough_v"].values)
        #         sd["slow_trough_norm_time"] = _slow_trough_norm_t(
        #             sd["threshold_t"].values,
        #             sd["slow_trough_t"].values,
        #             sd["trough_t"].values)

    def spikes_data(self):
        """ Return a list of spike feature dataframes"""
        return self._spikes_set[0]

    def sweep_features(self):
        return self._sweep_features.to_dict(orient='records')[0]

    def spike_ratio_features(self, features):
        dict_list = self._spikes_set[0].to_dict(orient='records')
        record = {}
        suffix = "_ratio"
        nspikes = len(dict_list)
        if nspikes==1:
            return
        for feature in features:
            last = dict_list[-1].get(feature)
            if last is None and nspikes >= 3:
                last = dict_list[-2].get(feature)
            if last is not None:
                record.update({feature+suffix: last/dict_list[0].get(feature)})

        return record