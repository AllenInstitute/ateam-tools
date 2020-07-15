#!/usr/bin/env python2
import ateam.sim.singlecell as sc
import os
import os.path
import numpy as np
from bmtk.simulator.bionet.biophys_params import BiophysParams
from ateam.sim.run.pbstools import run_hpc_bmtk

# import ateam.data.allensdk_tools as asdk
# cells_df = asdk.get_cells_ephys_df(manifest_file='/local1/storage/celltypes_data/manifest.json')
# cells_df = pd.read_csv('/home/tom.chartrand/projects/data/human_mouse_ephys_all_0127.csv', index_col=0)

networks_path = "/home/tom.chartrand/network/bmtk_networks/"
config_template = networks_path + "biophys_components_shared/default_config.json"
dynamics_path_base = networks_path + "biophys_components_shared/biophysical_neuron_templates/optim_param_{}.json"

# TODO: clean up old PARAMS_PATH issues
# params_dir = sc.OPT_PARAMS_PATH
sim_folder_path = networks_path + "singlecell/ih_psp_met/{}"
# sim_folder_path = networks_path + "singlecell/ih_psp_met_clamped/{}"

inh=False
n_duplicates=5
overwrite=True
sim_time=300
transient_time=500
edge_dict = {'syn_weight': 1e-4}

ih_factors = np.array([0.1, 0.5, 1, 2, 3])
props = ['genome.{}.gbar_Ih'.format(sec) for sec in ['soma', 'apic', 'dend']]

cells_list = [('frem', 774420848), ('frem', 786497938), 
    ('glp2r', 770268817), ('glp2r', 769228370)]

for ttype, cell_id in cells_list:
    dynamics_path = dynamics_path_base.format(cell_id)
    dp = BiophysParams.from_json(dynamics_path)
    
    ih_params = {prop: ih_factors*float(dp.get_nested(prop)) for prop in props}
    ih_params['ih_factor'] = ih_factors
    
    path = sim_folder_path.format(cell_id)
    sm = sc.build_epsp_batch(cell_id, path, overwrite=overwrite,
                          inh=inh, n_duplicates=n_duplicates, transient_time=transient_time,
                          edge_dict=edge_dict, node_dict={}, linked_dict=ih_params)
# add "--dryrun" to test
    run_hpc_bmtk( [sm.config_path, "-t", "00:40:00", "-p", "50", "--jobname", "hcem_singlecell_sim", "--conda", "bmtk"])
        