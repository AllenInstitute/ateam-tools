import numpy as np
import bmtk.analyzer.cell_vars as cell_vars
import matplotlib.pyplot as plt

def plot_v(config_file, ax=None, colors=None, gids=None, stack_sep=40, show=False):
    plot_var('v', config_file, ax=ax, colors=colors, gids=gids, stack_sep=stack_sep, show=show)

def plot_var(var_name, config_file, ax=None, colors=None, gids=None, stack_sep=None, compartments='origin', show=False):
    ax = ax or plt.gca()
    var_report = get_cellvar_report(config_file)
    tt = var_report.time_trace
    gids = gids or np.array(var_report.gids)

    for i, gid in enumerate(gids):
        color = colors[i] if colors else 'b'
        if stack_sep:
            ax.plot(tt, var_report.data(gid, var_name, compartments=compartments) + i*stack_sep, color=color)
        else:
            ax.plot(tt, var_report.data(gid, var_name, compartments=compartments), color=color)
    if show:
        plt.show()
                
def get_cellvar_report(config_file):
    report_file = cell_vars.get_cell_report(config_file)
    var_report = cell_vars.CellVarsFile(report_file)
    return var_report

def data_multicell(var_report, gid_list, var_name='v', time_window=None, compartments='origin'):
    data_iter = [var_report.data(gid, var_name, time_window, compartments) for gid in gid_list]
    data_stacked = np.stack(data_iter)
    return data_stacked