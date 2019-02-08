"""Functions for testing the synaptic input response of single biophysically detailed models.
Runs and analyzes batch simulations for two scenarios: single EPSPs and spiking response to Poisson input.
"""
import ateam.sim.setup.batch_builder as bb
import ateam.sim.setup.default_props as defaults
import ateam.sim.setup as sim
import ateam.sim.run.run_hpc as run_hpc
import numpy as np
import os

CONFIG_TEMPLATE_PATH = "/allen/aibs/mat/tmchartrand/bmtk_networks/biophys_components_shared/default_config.json"

def build_epsp_batch(cell_id, sim_folder, cell_name=None, inh=False, dmax=350, n_duplicates=10, edge_dict={}, node_dict={}):
    edge_props_base = {
        'nsyns': 5,
        'syn_weight': 1e-4,
        'delay': 0,
        'dynamics_params': 'AMPA_ExcToInh.json' if inh else 'AMPA_ExcToExc.json',
        'model_template': 'exp2syn'
    }
    edge_props_base.update(edge_dict)
    node_props = {
            'cell_name': cell_name or cell_id,
            'morphology': '{}.swc'.format(cell_id),
            'dynamics_params': 'optim_param_{}.json'.format(cell_id),
            'model_type': 'biophysical',
            'model_template': 'ctdb:Biophys1.hoc',
            'model_processing': 'aibs_allactive'
        }
    node_props.update(node_dict)
    config_path = os.path.join(sim_folder, "config.json")
    template = CONFIG_TEMPLATE_PATH
    sm = sim.SimManager.from_template(config_template=template, overwrite=True, config_path=config_path)
    opt_params_path = "/allen/aibs/mat/ateam_shared/All_active_params"
    sm.config.update_nested(components={"biophysical_neuron_models_dir": opt_params_path})

    input_net = bb.build_input_net_simple()

    distance_list = np.linspace(0, dmax, 10)
    interval = 5
    spike_time = 0.3
    linked_edge_props = {'distance_range_min':distance_list, 'distance_range_max':distance_list+interval}
    ind_edge_props = {'target_sections': ['a','b']}

    net = bb.build_batch_edge_props(input_net, node_props, edge_props_base, n_duplicates=n_duplicates,
                                    linked_dicts=[linked_edge_props], **ind_edge_props)

    sm.add_networks([net, input_net])
    sm.add_membrane_report()
    sm.write_spikeinput_vector(input_net.name, [spike_time])
    sm.save_network_files()
    return sm

def build_rates_batch(cell_id, sim_folder, cell_name=None, inh=False, max_input=2000, n_duplicates=10, edge_dict={}, node_dict={}, num_input=1):
    edge_props = {
        'nsyns': 50,
        'syn_weight': 1e-5, 
        'distance_range_min': 0,
        'distance_range_max': 200,
        'delay': 0,
        'target_sections': ['d','s'] if inh else ['a','b','s'],
        'dynamics_params': 'AMPA_ExcToInh.json' if inh else 'AMPA_ExcToExc.json',
        'model_template': 'exp2syn'
    }
    edge_props.update(edge_dict)
    config_path = os.path.join(sim_folder, "config.json")
    template = CONFIG_TEMPLATE_PATH
    sm = sim.SimManager.from_template(config_template=template, overwrite=True, config_path=config_path)
    opt_params_path = "/allen/aibs/mat/ateam_shared/All_active_params"
    sm.config.update_nested(components={"biophysical_neuron_models_dir": opt_params_path})

    node_props = {
            'cell_name': cell_name or cell_id,
            'morphology': '{}.swc'.format(cell_id),
            'dynamics_params': 'optim_param_{}.json'.format(cell_id),
            'model_type': 'biophysical',
            'model_template': 'ctdb:Biophys1.hoc',
            'model_processing': 'aibs_allactive'
        }
    node_props.update(node_dict)
    input_props = {'input_rate': np.linspace(1, max_input, 10), 'num_input': num_input}

    sm.add_membrane_report()
    net = bb.build_batch_all(sm, node_props, edge_props, input_props, n_duplicates=n_duplicates)
    return sm

def copy_sim_for_new_cell(sm, cell_id, folder, overwrite=False):
    """Save a copy of a single-cell simulation and update the cell ID in the config.
    
    Arguments:
        sm {SimManager} -- singlecell simulation to copy (must have saved network files)
        cell_id {int or str} -- Cell specimen ID for morph and params files
        folder {str} -- path to new file
    
    Keyword Arguments:
        overwrite {bool} -- (default: {False})
    """
    sm_mod = sm.save_copy(folder, overwrite)
    if sm_mod:
        node_props = {
                    'cell_name': cell_id,
                    'morphology': '{}.swc'.format(cell_id),
                    'dynamics_params': 'optim_param_{}.json'.format(cell_id)
        }
        sm_mod.update_node_type_props("batch", node_props)

def build_singlecell_sims(cells_list, base_path, inh=False):
    sim_path = os.path.join(base_path, "{cell}", "{sim}")
    cells_list = list(cells_list)
    base_id = cells_list.pop(0)
    sm_rates = build_rates_batch(base_id, sim_path.format(cell=base_id, sim='rate'),
                                 inh=inh, n_duplicates=10)
    sm_psp = build_epsp_batch(base_id, sim_path.format(cell=base_id, sim='psp'),
                              inh=inh, n_duplicates=10)
    for cell_id in cells_list:
        if cell_id != base_id:
            copy_sim_for_new_cell(sm_rates, cell_id, sim_path.format(cell=cell_id, sim='rate'), overwrite=False)
            copy_sim_for_new_cell(sm_psp, cell_id, sim_path.format(cell=cell_id, sim='psp'), overwrite=False)
            
def submit_singlecell_sims(cells_list, base_path, rerun=False):
    sim_path = os.path.join(base_path, "{cell}", "{sim}")
    for cell_id in cells_list:
        args_list = [os.path.join(sim_path.format(cell=cell_id, sim='rate'), "config.json"),
                    "-n2", "-t", "00:15:00"]
        if not rerun: args_list.append("--norerun")
        run_hpc.main(args_list)
        args_list = [os.path.join(sim_path.format(cell=cell_id, sim='psp'), "config.json")]
        if not rerun: args_list.append("--norerun")
        run_hpc.main(args_list)

import matplotlib.pyplot as plt
import pandas as pd
import os.path
import os
import ateam.analysis.nodes as nodes
import ateam.analysis.batch_analysis as ba
import ateam.analysis.cell_vars.psp_analysis as psp
import ateam.analysis.spikes as spikes

def extract_rates_data(config_path):
    sm = sim.SimManager(config_path)
    net = 'batch'
    nodes_df = nodes.create_node_table(sm.nodes_file(net), sm.node_types_file(net))
    df = nodes_df.join(pd.Series(spikes.get_rates_config(sm.config_path), name='rate'))
    df.fillna(0, inplace=True)
    return df

def extract_psp_data(config_path):
    sm = sim.SimManager(config_path)
    net = 'batch'
    nodes_df = nodes.create_node_table(sm.nodes_file(net), sm.node_types_file(net))
    spike_time = 0.3
    props = psp.epsp_analysis(config_file=sm.config_path, t_min=100, t_duration=1000, t_spike=1000*spike_time)
    syn_df = pd.DataFrame(props)
    syn_df.set_index('gid', inplace=True)
    df = syn_df.join(nodes_df)
    return df

def plot_rates_config(config_path, ax=None):
    df = extract_rates_data(config_path)
    df_stats = ba.plot_df_stats_comparison(df, 'input_rate', 'rate', 'target_sections', ax=ax)
    plt.legend(['apical','basal', 'somatic'])#, title='Target sections'
#     plt.gca().get_legend().remove()
    plt.xlabel('Input strength (a.u.)')
    plt.ylabel('Firing rate (Hz)')
    
def plot_psp_df(cellname, df, propname, label=None, ax=None):
    ba.plot_df_stats_comparison(df,'distance_range_max', propname, 'target_sections', ax=ax)
#     plt.legend(['apical','basal'], title='Target sections')
    ax.get_legend().remove()
    plt.xlabel('Distance from soma ($\mu$m)')
    label = label or propname
    plt.ylabel(label)

def plot_singlecell_rate_psp(cell, base_path=None, rate_sim=None, psp_sim=None, save_path=None, show=False):
    rate_sim = rate_sim or os.path.join(base_path, cell, "rate", "config.json")
    psp_sim = psp_sim or os.path.join(base_path, cell, "psp", "config.json")
    fig, axes = plt.subplots(2,2)
    plot_rates_config(rate_sim, ax=axes[0,0])

    df = extract_psp_data(psp_sim)
    plot_psp_df(cell, df, 'amp', ax=axes[0,1])
    plot_psp_df(cell, df, 'area', ax=axes[1,0])
    plot_psp_df(cell, df, 'delay', ax=axes[1,1])
    fig.suptitle("Cell {}".format(cell))
    plt.tight_layout()
    
    if save_path:
        save_dir = os.path.dirname(save_path)
        if not os.path.isdir(save_dir):
            os.mkdir(save_dir)
        plt.savefig(save_path)
    return fig

def save_plots_all_cells(cells, base_path, filename="all_cells.pdf"):
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
                plt.close()