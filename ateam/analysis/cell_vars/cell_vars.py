import numpy as np
import bmtk.analyzer.cell_vars as cell_vars
import matplotlib.pyplot as plt

def plot_v(config_file, ax=None, colors=None, gids=None, stack_sep=40, show=False):
    plot_var('v', config_file, ax=None, colors=None, gids=None, stack_sep=stack_sep, show=False)

def plot_var(var_name, config_file, ax=None, colors=None, gids=None, stack_sep=None, show=False):
    ax = ax or plt.gca()
    report_file = cell_vars.get_cell_report(config_file)
    var_report = cell_vars.CellVarsFile(report_file)
    tt = var_report.time_trace
    gids = gids or np.array(var_report.gids)

    for i, gid in enumerate(gids):
        color = colors[i] if colors else 'b'
        if stack_sep:
            ax.plot(tt, var_report.data(gid, var_name) + i*stack_sep, color=color)
        else:
            ax.plot(tt, var_report.data(gid, var_name), color=color)
    if show:
        plt.show()
                

def data_multicell(var_report, gid_list, var_name, time_window=None, compartments='origin'):
    data_iter = (var_report.data(gid, var_name, time_window, compartments) for gid in gid_list)
    data_stacked = np.stack(data_iter)
    return data_stacked