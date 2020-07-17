import numpy as np
import pandas as pd
import scipy.optimize as opt
import bmtk.analyzer.cell_vars as cell_vars
from .cell_vars import data_multicell

def remove_transient(var_report, t_spike, t_pre=10, t_max=1000):
    t0 = t_spike - t_pre
    v_all = data_multicell(var_report, gids, 'v', time_window=[t0, t_max])
    tt = np.arange(t0, var_report.t_stop, var_report.dt) - t_spike
    t_transient = tt[tt<0]
    f_exp = lambda t, a, b, tau: a*np.exp(-t/tau) + b
    jac_exp = lambda t, a, b, tau: np.stack([np.exp(-t/tau), np.ones_like(t), (a*t/tau**2)*np.exp(-t/tau)], axis=-1)
    for i in range(len(gids)):
        v_transient = v_all[i, tt<0]
        p0 = [0.1, v_transient[-1], 500]
        popt, _ = opt.curve_fit(f_exp, t_transient, v_transient, p0=p0, jac=jac_exp)
        v_all[i,:] -= f_exp(tt, *popt)
    return v_all, tt

def epsp_analysis(report_file=None, config_file='config.json', t_spike=None, t_max=1000, remove_transient=False):
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
    # TODO: check on t data in var report (start/stop off?)
    gids = np.array(var_report.gids)

    if remove_transient:
        v_all = v_all[:,tt>0]
        tt = tt[tt>0]
    else:     
        t0 = t_spike
        v_all = data_multicell(var_report, gids, 'v', time_window=[t0, t_max])
        # v array has shape (cells, time)
        tt = np.arange(t0, var_report.t_stop, var_report.dt) - t_spike
        v_all -= v_all[:,0,np.newaxis]

    i_peak = np.argmax(v_all, axis=-1)
    v_peak = np.max(v_all, axis=-1)
    cells_epsp = np.flatnonzero((i_peak != (v_all.shape[1]-1)) & (i_peak != 0) &
            (amp_min < v_peak) & (v_peak < amp_max))

    n_conn = len(cells_epsp)
    i_50_0 = np.zeros(n_conn, dtype=int)
    i_50_1 = np.zeros(n_conn, dtype=int)
    i_01 = np.zeros(n_conn, dtype=int)
    i_20 = np.zeros(n_conn, dtype=int)
    i_80 = np.zeros(n_conn, dtype=int)
    area = np.zeros(n_conn)
    for i, cell in enumerate(cells_epsp):
        i_up = np.flatnonzero(v_all[cell,:] > 0.5*v_peak[cell])
        i_50_0[i] = i_up[0]
        i_50_1[i] = i_up[-1]
        i_01[i] = np.flatnonzero(v_all[cell,:] > 0.01*v_peak[cell])[0]
        i_20[i] = np.flatnonzero(v_all[cell,:] > 0.2*v_peak[cell])[0]
        i_80[i] = np.flatnonzero(v_all[cell,:] > 0.8*v_peak[cell])[0]

        area[i] = np.sum(v_all[cell, :])*var_report.dt

    t_start = tt[i_01]
    peak_time = tt[i_peak[cells_epsp]] - t_start
    return {
        'gid': gids[cells_epsp], 
        'amp': v_peak[cells_epsp], 
        'width': tt[i_50_1] - tt[i_50_0], 
        'latency': t_start, 
        'peak_time': peak_time,
        'rise_time': tt[i_80] - tt[i_20], 
        'area': area,
        }
        
