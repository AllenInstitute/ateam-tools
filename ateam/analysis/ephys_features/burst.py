import numpy as np
import warnings
import ateam.data.lims as lims
from ateam.data.lims import LimsReader
from allensdk.core.nwb_data_set import NwbDataSet
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor
from ateam.analysis.ephys_features.triblip import postspike_v_all
from ateam.analysis.ephys_features.common import stim_props

def cost_mean(x):
    mu = np.mean(x)
    ss = np.sum((x - mu) ** 2)
    return ss


def changepoint_detect(input_data, penalty=10., cost_func=cost_mean):
    data = input_data.copy()
    if len(data.shape) == 1: # convert to 2D equivalent if vector
        data = data.reshape((len(data), -1))

    n = data.shape[0]

    checklist = np.zeros(n + 1, dtype=int)
    lastchange_likelihood = np.zeros(n + 1, dtype=float)
    lastchange_cpts = np.zeros(n + 1, dtype=int)
    n_cpts = np.zeros(n + 1, dtype=int)
    temp_likelihood = np.zeros(n + 1, dtype=float)
    temp_t = np.zeros(n + 1, dtype=float)

    lastchange_likelihood[0] = -penalty
    lastchange_cpts[0] = 0
    n_cpts[0] = 0

    lastchange_likelihood[1] = cost_func(data[:1, :])
    lastchange_cpts[1] = 0
    n_cpts[1] = 1

    n_checklist = 2
    checklist[0] = 0
    checklist[1] = 1

    for tstar in range(2, n + 1):
        if lastchange_likelihood[tstar] == 0:
            for i in range(n_checklist):
                temp_likelihood[i] = (lastchange_likelihood[checklist[i]] +
                                      cost_func(data[checklist[i]:tstar, :]) +
                                      penalty)
            which_out = np.argmin(temp_likelihood[:n_checklist])
            min_out = temp_likelihood[which_out]
            lastchange_likelihood[tstar] = min_out
            lastchange_cpts[tstar] = checklist[which_out]
            n_cpts[tstar] = n_cpts[lastchange_cpts[tstar]] + 1

            n_check_temp = 0
            for i in range(n_checklist):
                if temp_likelihood[i] <= lastchange_likelihood[tstar] + penalty:
                    checklist[n_check_temp] = checklist[i]
                    n_check_temp += 1
            n_checklist = n_check_temp

        checklist[n_checklist] = tstar
        n_checklist += 1

    last = n
    cpts_out = []
    while last != 0:
        cpts_out.insert(0, last)
        last = lastchange_cpts[last]
    return cpts_out[:-1]


def process_cell(specimen_id, passed_only=True, use_silence=True):
    lr = LimsReader()
    nwb_path = lr.get_nwb_path_from_lims(specimen_id)
    longs = lr.get_sweeps(specimen_id, "Long Square")

    dataset = NwbDataSet(nwb_path)
    results = []
    segments_list = []
    try:
        for l in longs:
            v, i, t = lims.get_sweep_v_i_t_from_set(dataset, l, window_data=False)
            stim, amp = stim_props(v, i, t)
            
            if amp < 0:
                continue
            sweepex = EphysSweepFeatureExtractor(t=t, v=v, i=i, start=stim[0], end=stim[1])  
            sweepex.process_spikes()
            st = sweepex.spike_feature("peak_t")

            if len(st) < 3:
                segments_list.append(st)
                continue

            segment_info = segment(st)
            segments = segment_info["segments"]

            n_values = [s["n"] for s in segments]
            rate_values = [s["avg_rate"] for s in segments]


            n_max_ind = np.argmax(n_values)
            base_rate = rate_values[n_max_ind]
            max_rate = np.max(rate_values)
            burst_ratio = max_rate / base_rate
            # pause_ratio = base_rate / np.min(rate_values)
            var_ratio = max_rate / np.min(rate_values)

            burst_index_n = 1 - 1./burst_ratio

            t_max_ind = np.argmax([s["time"] for s in segments])
            burst_ratio_t = max_rate / rate_values[t_max_ind]
            burst_index_t = 1 - 1./burst_ratio_t

            t_silent = stim[1] - st[-1]
            burst_index_init = burst_index_t
            if t_silent > segments[t_max_ind]["time"]:
                burst_index_init = 1

            # if use_silence:
            #     burst_index_t = burst_index_init 

                # segments.append({
                #     "time": t_silent,
                #     "n": 0,
                #     "avg_rate": 0,
                # })

            vpost = postspike_v_all(sweepex, delta_t=0.005)
            vpost_first = vpost[ segments[0]["n"] - 1]

            vahp = sweepex.spike_feature("trough_v")
            vahp_first = np.mean(vahp) if len(segments)==1 else vahp[ segments[0]["n"] - 1]

            sweep_dict = {
                "cell_id": specimen_id,
                "sweep_num": l,
                "n_burst": n_values[np.argmax(rate_values)],
                # "stimulus_amplitude": amp,
                "burst_ratio": burst_ratio,
                "burst_ratio_t": burst_ratio_t,
                "burst_ratio_var": var_ratio,
                # "pause_ratio": pause_ratio,
                "burst_index_init": burst_index_init,
                "burst_index_n": burst_index_n,
                "burst_index_t": burst_index_t,

                "vpost": np.mean(vpost),
                "delta_vpost": vpost_first - np.mean(vpost),

                "vahp": np.mean(vpost),
                "delta_vahp": vahp_first - np.mean(vahp)
            }
            results.append(sweep_dict)
            segments_list.append(segments)
    except ValueError as e:
        print("Exception on", specimen_id)
        print(e)
        # raise

    return results


def segment(spike_times, cv_threshold=0.18, stim_end_time=None):
    penalty_set = (10 ** np.linspace(0., 4., 50)) / 1e5
    isis = np.diff(spike_times)
    inst_freq = 1. / isis

    for penalty in reversed(penalty_set):
        cpts = changepoint_detect(isis, penalty=penalty)

        segments_all_pass = True
        for start, end in zip([0] + cpts, cpts + [len(isis)]):
            if end - start == 0:
                continue
            seg_freq = inst_freq[start:end]
            cv = seg_freq.std() / seg_freq.mean()
            if cv > cv_threshold:
                segments_all_pass = False
                break
        if segments_all_pass:
            break

    segments = []
    for start, end in zip([0] + cpts, cpts + [len(isis)]):
        seg_freq = inst_freq[start:end]
        segments.append({
            "start": start,
            "end": end,
            "time": spike_times[end] - spike_times[start],
            "n": end - start,
            "avg_rate": seg_freq.mean(),
        })
    

    return {"penalty": penalty, "segments": segments}