#!/usr/bin/env python2
import os.path
import numpy as np
import pandas as pd
import ateam.sim.setup.batch_builder as bb
from ateam.sim.run.pbstools import run_hpc_bmtk
from bmtk.simulator.bionet.biophys_params import BiophysParams
import ateam.sim.singlecell as sc
import ateam.sim.setup.default_props as defaults
from ateam.sim.setup import SimManager

# import ateam.data.allensdk_tools as asdk
# cells_df = asdk.get_cells_ephys_df(manifest_file='/local1/storage/celltypes_data/manifest.json')
# rheo = cells_df.
cells_df = pd.read_csv('/home/tom.chartrand/projects/data/human_mouse_ephys_all_0127.csv', index_col=0)
rheo = cells_df.rheobase_i

networks_path = "/home/tom.chartrand/network/bmtk_networks/"
config_template = networks_path + "biophys_components_shared/default_config.json"
dynamics_path_base = networks_path + "biophys_components_shared/biophysical_neuron_templates/optim_param_{}.json"

# sim_path = networks_path + "singlecell/rheo_new_clamped/{}/config.json"

sim_path = networks_path + "singlecell/rheo_new/{}/config.json"
g_factors = np.array([0.1, 0.5, 1, 2, 3])

# sim_path = networks_path + "singlecell/rheo_new_fine/{}/config.json"
# g_factors = np.array([0.8, 0.9, 1, 1.2, 1.4])

g_props = {
    'cah': ['genome.{}.gbar_Ca_HVA'.format(sec) for sec in ['soma', 'axon']],
    'cal': ['genome.{}.gbar_Ca_LVA'.format(sec) for sec in ['soma', 'axon']],
    'kt': ['genome.{}.gbar_K_Tst'.format(sec) for sec in ['soma', 'axon']],
    'sk': ['genome.{}.gbar_SK'.format(sec) for sec in ['soma', 'axon']],
    'ih': ['genome.{}.gbar_Ih'.format(sec) for sec in ['soma', 'apic', 'dend']],
    'kv3': ['genome.{}.gbar_Kv3_1'.format(sec) for sec in ['soma', 'apic', 'dend', 'axon']],
}


# rheo_offsets = [-10, -5, 0, 5, 10, 50]
rheo_offsets = [0, 50]
sim_time = 1600.

cells_list = [('frem', 774420848), ('frem', 786497938), 
    ('glp2r', 770268817), ('glp2r', 769228370)]

for ttype, cell_id in cells_list:
    dynamics_path = dynamics_path_base.format(cell_id)
    dp = BiophysParams.from_json(dynamics_path)

    input_dict = {
        "amp": (rheo.loc[cell_id] + rheo_offsets)/1000.,
        "delay": 500.,
        "duration": 1000.
    }
    
    node_dict = defaults.cellprops_active(cell_id)
    # node_dict.update({'dynamics_params': dynamics_path})
    

    for channel in ['cal']:
    # for channel in ['ih', 'kt']:
        cell_name = "{}_{}".format( cell_id, channel)
        node_dict.update({'cell_name': cell_name})
        g_params = {prop: g_factors*float(dp.get_nested(prop)) for prop in g_props[channel]}
        g_params['g_factor'] = g_factors
        sm = SimManager.from_template(config_path=sim_path.format(cell_name), config_template=config_template, overwrite=True)
        net = bb.build_batch_all(sm, node_dict, {}, input_dict, linked_dicts=[g_params])
    
        sm.add_membrane_report()
        sm.sim_time = sim_time
        sm.save_network_files()

        # add "--dryrun" to test
        run_hpc_bmtk( [sm.config_path, "-t", "00:40:00", "-p", "20", "--jobname", "hcem_singlecell_sim", "--conda", "bmtk"])
        
        # need to use the direct one to debug
        # from ateam.sim.run import run_bionet
        # run_bionet.run(sm.config_path)