"""Functions for testing the synaptic input response of single biophysically detailed models.
Runs and analyzes batch simulations for two scenarios: single EPSPs and spiking response to Poisson input.
"""
import ateam.sim.setup.batch_builder as bb
import ateam.sim.setup.default_props as defaults
import ateam.sim.setup as sim
import numpy as np
import os

CONFIG_TEMPLATE_PATH = "/allen/aibs/mat/tmchartrand/bmtk_networks/biophys_components_shared/default_config.json"
opt_params_path = "/allen/aibs/mat/ateam_shared/All_active_params"

def get_hof_params_path(cell_id, hof_id):
    path = "/allen/aibs/mat/ateam_shared/Human_Model_Fit_Metrics/{cell_id}/fitted_params/hof_param_{cell_id}_{hof_id}.json"
    return path.format(cell_id=cell_id, hof_id=hof_id)

def edge_props_shared():
    return {
        'nsyns': 5,
        'syn_weight': 5e-5,
        'delay': 0,
        'dynamics_params': 'AMPA_bbp.json',
        'model_template': 'exp2syn'
        }

def get_node_props(cell_id, cell_name=None):
    node_props = {
            'cell_name': cell_name or cell_id,
            'cell_id': cell_id,
            'morphology': '{}.swc'.format(cell_id),
            'dynamics_params': 'optim_param_{}.json'.format(cell_id),
            'model_type': 'biophysical',
            'model_template': 'ctdb:Biophys1.hoc',
            'model_processing': 'aibs_allactive_ani'
        }
    return node_props

def build_epsp_batch(cell_id, sim_folder, cell_name=None, inh=False,
    dmax=400, n_duplicates=10, edge_dict={}, node_dict={}, overwrite=False):
    synapse_cluster_scale = 20
    sim_time = 300
    input_props = {'spike_time':0.2}
    edge_props = edge_props_shared()
    edge_props.update({
        # These must be included here if overwritten in linked props
        'target_sections': [0],
        'distance_range_min': [0],
        'distance_range_max': [0]
    })
    ndist = 10
    nsecs = 1 if inh else 2
    distance_min = np.linspace(0, dmax, ndist)
    distance_max = distance_min + synapse_cluster_scale
    secs = ['s'] + ndist*['d'] if inh else ['s'] + ndist*['a'] + ndist*['d']
    linked_edge_props = {'target_sections': secs,
        'distance_range_min':[0]+nsecs*list(distance_min), 
        'distance_range_max':[1]+nsecs*list(distance_max)}

    edge_props.update(edge_dict)
    node_props = get_node_props(cell_id, cell_name)
    node_props.update(node_dict)
    config_path = os.path.join(sim_folder, "config.json")
    sm = sim.SimManager.from_template(config_template=CONFIG_TEMPLATE_PATH, overwrite=True, config_path=config_path)
    sm.config.update_nested(components={"biophysical_neuron_models_dir": opt_params_path})
    sm.add_membrane_report()
    sm.sim_time = sim_time

    net = bb.build_batch_all(sm, node_props, edge_props, input_props, linked_dicts=[linked_edge_props], n_duplicates=n_duplicates, use_abs_paths=True)
    return sm

def build_rates_batch(cell_id, sim_folder, cell_name=None, inh=False,
    max_input=800, num_input=10, n_duplicates=10, edge_dict={}, node_dict={}, overwrite=False):
    sections = ['s','d'] if inh else ['s','a','d']
    
    edge_props = edge_props_shared()
    edge_props.update({
        'distance_range_min': 0,
        'distance_range_max': 1000,
        'target_sections': sections
    })
    nrates = 10
    input_props = {'input_rate': np.linspace(1, max_input, nrates), 'num_input': num_input}
    
    edge_props.update(edge_dict)
    config_path = os.path.join(sim_folder, "config.json")
    template = CONFIG_TEMPLATE_PATH
    sm = sim.SimManager.from_template(config_template=template, overwrite=True, config_path=config_path)
    sm.config.update_nested(components={"biophysical_neuron_models_dir": opt_params_path})

    node_props = get_node_props(cell_id, cell_name)
    node_props.update(node_dict)

    # sm.add_membrane_report()
    net = bb.build_batch_all(sm, node_props, edge_props, input_props, n_duplicates=n_duplicates, use_abs_paths=True)
    return sm

def build_hof_sims(base_path, hof_path, cell_id, num_sims=40, overwrite=False, level='folder'):
    def build_hof_internal(sim_type):
        sm = sim.SimManager(config_path(base_path, cell_id, sim_type))
        for hof_id in range(num_sims):
            folder = sim_path(hof_path, cell_id, sim_type, hof_id)
            copy_sim_for_hof(sm, cell_id, hof_id, folder, overwrite=overwrite, level=level)
    build_hof_internal("rate")
    build_hof_internal("psp")

def copy_sim_for_hof(sm, cell_id, hof_id, folder, overwrite=False, level='folder'):
    sm_mod = sm.save_copy(folder, overwrite)
    if sm_mod:
        node_props = {
                    'cell_id': cell_id,
                    'cell_name': "{}_hof_{}".format(cell_id, hof_id),
                    'morphology': '{}.swc'.format(cell_id),
                    'dynamics_params': get_hof_params_path(cell_id, hof_id)
        }
        sm_mod.update_node_type_props("batch", node_props)

def copy_sim_for_new_cell(sm, cell_id, folder, overwrite=False, level='folder'):
    """Save a copy of a single-cell simulation and update the cell ID in the config.
    
    Arguments:
        sm {SimManager} -- singlecell simulation to copy (must have saved network files)
        cell_id {int or str} -- Cell specimen ID for morph and params files
        folder {str} -- path to new file
    
    Keyword Arguments:
        overwrite {bool} -- (default: {False})
    """
    sm_mod = sm.save_copy(folder, overwrite=overwrite, level=level)
    if sm_mod:
        node_props = {
                    'cell_name': cell_id,
                    'morphology': '{}.swc'.format(cell_id),
                    'dynamics_params': 'optim_param_{}.json'.format(cell_id)
        }
        sm_mod.update_node_type_props("batch", node_props)

def build_singlecell_sims(cells_list, base_path, inh=False, n_duplicates=10, overwrite=False, level='folder'):
    sim_path = os.path.join(base_path, "{cell}", "{sim}")
    base_id = cells_list[0]
    sm_rates = build_rates_batch(base_id, sim_path.format(cell="base", sim='rate'),
                                 inh=inh, n_duplicates=n_duplicates, overwrite=overwrite)
    sm_psp = build_epsp_batch(base_id, sim_path.format(cell="base", sim='psp'),
                              inh=inh, n_duplicates=n_duplicates, overwrite=overwrite)
    for cell_id in cells_list:
        copy_sim_for_new_cell(sm_rates, cell_id, sim_path.format(cell=cell_id, sim='rate'), overwrite=overwrite, level=level)
        copy_sim_for_new_cell(sm_psp, cell_id, sim_path.format(cell=cell_id, sim='psp'), overwrite=overwrite, level=level)
            
def sim_path(base_path, cell, sim, hof_id=""):
    sim_path = os.path.join(base_path, str(cell), str(hof_id), sim)
    return sim_path

def config_path(base_path, cell, sim, hof_id=""):
    return os.path.join(sim_path(base_path, cell, sim, hof_id), "config.json")