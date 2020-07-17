import numpy as np
import bmtk.analyzer.cell_vars as cell_vars
import matplotlib.pyplot as plt
from bmtk.utils.cell_vars.var_reader import CellVarsFile

def plot_v(config_file, ax=None, colors=None, gids=None, stack_sep=40, time_window=None, show=False):
    plot_var('v', config_file, ax=ax, colors=colors, gids=gids, stack_sep=stack_sep, time_window=time_window, show=show)

def plot_var(var_name, config_file, ax=None, colors=None, gids=None, stack_sep=None, compartments='origin', time_window=None, show=False):
    ax = ax or plt.gca()
    var_report = get_cellvar_report(config_file)
    # tt = var_report.time_trace
    if gids is None:
        gids = np.array(var_report.gids)

    for i, gid in enumerate(gids):
        color = colors[i] if colors is not None else None
        v, tt = var_report.data_timeseries(gid, var_name, compartments=compartments, time_window=time_window)
        if stack_sep:
            ax.plot(tt, v + i*stack_sep, color=color)
        else:
            ax.plot(tt, v, color=color)
    if show:
        plt.show()
                
def get_cellvar_report(config_file):
    report_file = cell_vars.get_cell_report(config_file)
    var_report = CellVarsFile(report_file)
    return var_report

# TODO: move into bmtk cell_vars, could load more efficiently?
def data_multicell(var_report, gid_list=None, var_name='v', time_window=None, compartments='origin'):
    if gid_list is None:
        gid_list = var_report.gids
    data_iter = [var_report.data(gid, var_name, time_window, compartments) for gid in gid_list]
    data_stacked = np.stack(data_iter)
    return data_stacked