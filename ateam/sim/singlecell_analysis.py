
import matplotlib.pyplot as plt
import pandas as pd
import os.path
import os
from functools import partial
import ateam.analysis.bmtk.nodes as nodes
import ateam.analysis.bmtk.psp_analysis as psp
import ateam.analysis.bmtk.spikes as spikes
import ateam.analysis.dataframes as dfa
import ateam.sim.setup as sim
from ateam.sim.singlecell import sim_path, config_path

def fit_path(base_path, cells_list=None):
    cells_list = cells_list or os.listdir(base_path)
    fit_function = partial(fit_df_all, base_path=base_path)
    fit_all = map(fit_function, cells_list)
    cells_fit_df = pd.concat(fit_all, axis=0, keys=cells_list, names=['cell'])
    
    cells_fit_df = dfa.flatten_columns(cells_fit_df).reset_index(level="target_sections")
    return cells_fit_df

def fit_df_all(cell, base_path, extra=None, hof_id=""):
    # TODO: clean up inclusion of extra field
    psp = fit_df_psp(cell, base_path, extra=extra, hof_id=hof_id)
    rate = fit_df_rate(cell, base_path, hof_id=hof_id)
    if psp is None and rate is None:
        return None
    fit_df = pd.concat([psp, rate], axis=1, sort=True)
    fit_df.index.name = "target_sections"
    return fit_df

def fit_df_rate(cell, base_path, hof_id="", extra=None):
    try:
        df = extract_rates_data( config_path(base_path, cell, "rate", hof_id=hof_id) )
    except Exception as e:
        print('Caught exception in cell (x = {}):'.format(cell)) 
        # traceback.print_exc()
        print(e)
        return None
    
    compare = "target_sections"
    xvar = "input_rate"
    yvar = 'rate'
    # extra = ['morphology']
    functions = [dfa.linfit, dfa.threshold]
    grouped_df = df.groupby([compare])
    
    combined_fcn = dfa.combine_functions_hierarchical(xvar, yvar, functions, dropna=True)
    fit_df = grouped_df.apply(combined_fcn)
    
    fit_df = pd.concat([fit_df], axis=1, keys=[yvar], names=['prop'])
    if extra:
        fit_df = pd.concat([grouped_df[extra].agg(dfa.summary), fit_df], axis=1, sort=True)
    return fit_df

def fit_df_psp(cell, base_path, hof_id="", extra=None):
    try:
        df = extract_psp_data( config_path(base_path, cell, "psp", hof_id=hof_id) )
    except Exception as e:
        print('Caught exception in cell (x = {}):'.format(cell)) 
        # traceback.print_exc()
        print(e)
        return None
    
    compare = "target_sections"
    xvar = "distance_range_min"
    ylist = ['amp','delay','area']
    grouped_df = df.groupby([compare])
    
    functions = [dfa.linfit]
    fit_list = [grouped_df.apply(dfa.combine_functions_hierarchical(xvar, yvar, functions, dropna=True)) 
               for yvar in ylist]
    
    fit_df = pd.concat(fit_list, axis=1, keys=ylist, names=['prop'])
    if extra:
        fit_df = pd.concat([grouped_df[extra].agg(dfa.summary), fit_df], axis=1, sort=True)
    return fit_df

def extract_rates_data(config_path):
    sm = sim.SimManager(config_path)
    net = 'batch'
    nodes_df = nodes.create_node_table(sm.nodes_file(net), sm.node_types_file(net))
    df = nodes_df.join(pd.Series(spikes.get_rates_config(sm.config_path), name='rate'))
    return df

def extract_psp_data(config_path):
    sm = sim.SimManager(config_path)
    net = 'batch'
    nodes_df = nodes.create_node_table(sm.nodes_file(net), sm.node_types_file(net))
    spike_time = nodes_df.spike_time.unique()
    assert(len(spike_time)==1)
    spike_time = spike_time[0]
    # TODO: work with vector of spike times?
    props = psp.epsp_analysis(config_file=sm.config_path, t_min=100, t_duration=1000, t_spike=1000*spike_time)
    syn_df = pd.DataFrame(props)
    syn_df.set_index('gid', inplace=True)
    df = syn_df.join(nodes_df)
    return df

def plot_rates_config(config_path, ax=None):
    ax = ax or plt.axes()
    plt.sca(ax)
    df = extract_rates_data(config_path)
    df.fillna(0, inplace=True)
    df_stats = plot_df_stats_comparison(df, 'input_rate', 'rate', 'target_sections', ax=ax)
    # plt.legend(['apical','basal', 'somatic'])#, title='Target sections'
    # plt.gca().get_legend().remove()
    plt.xlabel('Input strength (a.u.)')
    plt.ylabel('Firing rate (Hz)')
    
def plot_psp_df(cellname, df, propname, label=None, ax=None, legend=False):
    ax = ax or plt.axes()
    plt.sca(ax)
    plot_df_stats_comparison(df,'distance_range_max', propname, 'target_sections', ax=ax)
    if not legend:
        ax.get_legend().remove()
    plt.xlabel('Distance from soma ($\mu$m)')
    label = label or propname
    plt.ylabel(label)
    plt.ylim(ymin=0)


def plot_singlecell_rate_psp(cell, base_path=None, rate_sim=None, psp_sim=None, save_path=None, show=False):
    rate_sim = rate_sim or os.path.join(base_path, cell, "rate", "config.json")
    psp_sim = psp_sim or os.path.join(base_path, cell, "psp", "config.json")
    fig, axes = plt.subplots(2,2)
    try:
        plot_rates_config(rate_sim, ax=axes[0,0])
    except:
        pass
    try:
        df = extract_psp_data(psp_sim)
        plot_psp_df(cell, df, 'amp', label='amp (mV)', ax=axes[0,1])
        plot_psp_df(cell, df, 'area', label='area (mV*ms)', ax=axes[1,0])
        plot_psp_df(cell, df, 'delay', label='delay (ms)', ax=axes[1,1])
    except:
        pass
    fig.suptitle("Cell {}".format(cell))
    plt.tight_layout()
    
    if save_path:
        save_dir = os.path.dirname(save_path)
        if not os.path.isdir(save_dir):
            os.mkdir(save_dir)
        plt.savefig(save_path)
    return fig

def save_plots_all_cells(base_path, cells=None, filename="all_cells.pdf"):
    # import traceback
    cells = cells or os.listdir(base_path)
    if not os.path.isabs(filename):
        filename = os.path.join(base_path, filename)
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(filename) as pdf:
        for cell in cells:
#             file_path = os.path.join(fig_folder, '{}.pdf'.format(cell))
            try:
                fig = plot_singlecell_rate_psp(cell, base_path)
                pdf.savefig(fig)
                plt.close()
            except Exception as exc:
                print("Cell {} failed.\n".format(cell), exc)
                # traceback.print_exc()
                plt.close()

# Note: this might all be in the seaborn functionality, and better documented!

def plot_df_stats_comparison_multiple(df, xvar, yvar, compare, separate, figsize=None):
    df_agg = df.groupby([xvar, compare, separate])[yvar].agg(['mean','std'])
    subsets = df_agg.index.unique(level=separate)
    nplot = len(subsets)
    fig, ax = plt.subplots(nplot, figsize=figsize)
    for i, name in enumerate(subsets):
        df_sub = df_agg.xs(name, level=separate)
        for name_compare in df_agg.index.unique(level=compare).sort_values():
            df_sub.xs(name_compare, level=compare).plot(y='mean', yerr='std', ax=ax[i], label=name_compare)
        ax[i].set(ylabel=yvar, title=name)
        ax[i].legend(title=compare)
    plt.tight_layout()

def plot_df_stats_comparison(df, xvar, yvar, compare, ax=None, **plot_args):
    df_agg = df.groupby([xvar, compare])[yvar].agg(['mean','std'])
    ax = ax or plt.axes()
    for name_compare in df_agg.index.unique(level=compare).sort_values():
        df_agg.xs(name_compare, level=compare).plot(y='mean', yerr='std', ax=ax, label=name_compare, **plot_args)
    ax.set(ylabel=yvar)
    ax.legend(title=compare)
    return df_agg