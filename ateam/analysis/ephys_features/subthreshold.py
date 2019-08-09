import numpy as np
import ateam.data.lims as lims
import allensdk.ephys.ephys_features as ft
from allensdk.core.nwb_data_set import NwbDataSet
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor
from ateam.analysis.ephys_features.common import stim_props

def process_cell(specimen_id, passed_only=True):
    lr = lims.LimsReader()
    nwb_path = lr.get_nwb_path_from_lims(specimen_id)
    longs = lr.get_sweeps(specimen_id, "Long Square")

    dataset = NwbDataSet(nwb_path)
    results = []
    # voltage threshold for quick spiking check
    thresh = 0 
    try:
        for l in longs:
            v, i, t = lims.get_sweep_v_i_t_from_set(dataset, l, window_data=False)
            stim, amp = stim_props(v, i, t)
            
            if amp < 0 or any(v > thresh):
                continue
            sweepex = EphysSweepFeatureExtractor(t=t, v=v, i=i, start=stim[0], end=stim[1])  
            # want to save absolute not only fraction, and modify deflection params
            # sag = sweepex.estimate_sag()
            # v_peak, peak_index = sweepex.voltage_deflection()
            start, end = stim
            peak_width=0.005
            peak_interval = 0.05
            start_index = ft.find_time_index(t, start)
            int_index = ft.find_time_index(t, start+peak_interval)
            peak_index = start_index + np.argmax(v[start_index:int_index])
            v_peak_avg = ft.average_voltage(v, t, start=t[peak_index] - peak_width / 2.,
                                            end=t[peak_index] + peak_width / 2.)
            v_baseline = sweepex.sweep_feature("v_baseline")
            v_steady = ft.average_voltage(v, t, start=end - sweepex.baseline_interval, end=end)
            hump_size = (v_peak_avg - v_steady)
            deflection = (v_peak_avg - v_baseline)
            hump_norm = hump_size / deflection

            sweep_dict = {
                "cell_id": specimen_id,
                "sweep_num": l,
                "current": amp,
                "hump": hump_size,
                "hump_norm": hump_norm,
                "deflection": deflection
            }
            results.append(sweep_dict)
    except ValueError as e:
        print("Exception on", specimen_id)
        print(e)
        # raise

    return results
    
