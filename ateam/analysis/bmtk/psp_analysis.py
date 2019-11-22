import numpy as np
import pandas as pd
import bmtk.analyzer.cell_vars as cell_vars
from .cell_vars import data_multicell

def epsp_analysis(report_file=None, config_file='config.json', t_spike=None, t_max=1000):
    """ Analyze epsp's as contained in a bionet cell vars report,
    for a simulation consisting of a single presynaptic spike,
     together with postsynaptic effects.
    Specify report_file directly, or lookup in simulation config file.

    If exact presynaptic input time is known, specify as t_spike
    Otherwise, include presynaptic cell and the method 
    will search for a spike starting at t_min.
    """
    amp_min = 0.05
    amp_max = 50
    report_file = report_file or cell_vars.get_cell_report(config_file)
    var_report = cell_vars.CellVarsFile(report_file)
    gids = np.array(var_report.gids)

    v_all = data_multicell(var_report, gids, 'v', time_window=[t_spike, t_max])
    # v array has shape (cells, time)
    
    tt = np.arange(t_spike, t_max, var_report.dt) - t_spike

    i_peak = np.argmax(v_all, axis=-1)
    v_peak = np.max(v_all, axis=-1)
    v_bl = v_all[:,0]
    amp = v_peak - v_bl
    cells_epsp = np.flatnonzero((i_peak != (v_all.shape[1]-1)) & (i_peak != 0) &
            (amp_min < amp) & (amp < amp_max))

    v_mid = (v_bl + v_peak)/2
    v_10 = v_bl + 0.1*(v_peak-v_bl)
    n_conn = len(cells_epsp)
    i_start = np.zeros(n_conn, dtype=int)
    i_end = np.zeros(n_conn, dtype=int)
    area = np.zeros(n_conn)
    for i, cell in enumerate(cells_epsp):
        i_up = np.flatnonzero(v_all[cell,:] > v_mid[cell])
        i_start[i] = i_up[0]
        i_end[i] = i_up[-1]

        i_area = np.flatnonzero(v_all[cell,:] > v_10[cell])
        area[i] = np.sum(v_all[cell, i_area] - v_bl[cell])*var_report.dt

    t_start = tt[i_start]
    t_end = tt[i_end]
    width = t_end - t_start
    return {
        'gid': gids[cells_epsp], 
        'amp': amp[cells_epsp], 
        'width': width, 
        'delay': t_start, 
        'peak_delay': tt[i_peak[cells_epsp]], 
        'area': area,
        }
        
