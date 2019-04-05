import numpy as np
import warnings
import ateam.data.lims as lims
from ateam.data.lims import LimsReader
from allensdk.core.nwb_data_set import NwbDataSet
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor

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

def process_cell(specimen_id, passed_only=True, use_silence=False):
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
            st = spike_times(v, i, t, start=stim[0], end=stim[1])

            if len(st) < 2:
                segments_list.append(st)
                continue

            segment_info = segment(st)
            segments = segment_info["segments"]

            if use_silence:
                # if t_silent > 2*isi_last:
                t_silent = stim[1] - st[-1]
                segments.append({
                    "time": t_silent,
                    "n": 0,
                    "avg_rate": 0,
                })

            n_values = [s["n"] for s in segments]
            rate_values = [s["avg_rate"] for s in segments]


            n_max_ind = np.argmax(n_values)
            base_rate = rate_values[n_max_ind]
            max_rate = np.max(rate_values)
            burst_ratio =  max_rate/ base_rate
            # pause_ratio = base_rate / np.min(rate_values)

            burst_index_n = 1 - base_rate / max_rate
            t_max_ind = np.argmax([s["time"] for s in segments])
            burst_index_t = 1 - rate_values[t_max_ind] / max_rate

            sweep_dict = {
                "cell_id": specimen_id,
                "sweep_num": l,
                # "stimulus_amplitude": amp,
                "burst_ratio": burst_ratio,
                # "pause_ratio": pause_ratio,
                "burst_index_n": burst_index_n,
                "burst_index_t": burst_index_t
            }
            results.append(sweep_dict)
            segments_list.append(segments)
    except:
        print "Exception on", specimen_id
        raise

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