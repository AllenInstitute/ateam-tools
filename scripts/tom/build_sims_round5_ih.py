import ateam.sim.singlecell as sc
import os
import os.path
import numpy as np
import ateam.data.allensdk_tools as asdk
from bmtk.simulator.bionet.biophys_params import BiophysParams

# April 17, 2020
# Changes: vary Ih with longer transient

# bootstrap from initial run
cells_path = "/allen/aibs/mat/tmchartrand/bmtk_networks/singlecell/round1"
cells_exc = os.listdir(os.path.join(cells_path,"PC"))
cells_list = [int(cell) for cell in cells_exc]

base_path = "/allen/aibs/mat/tmchartrand/bmtk_networks/singlecell/round5/Ih"
params_dir = sc.OPT_PARAMS_PATH

inh=False
n_duplicates=5
overwrite=True
sim_time=300
transient_time=500

sim_path = os.path.join(base_path, "{cell}", "{sim}")
base_id = cells_list[0]
edge_dict = {'syn_weight': 1e-4}

ih_factors = np.array([0.1, 0.5, 1, 2, 5])
props = ['genome.{}.gbar_Ih'.format(sec) for sec in ['soma', 'apic', 'dend']]

cells_df = asdk.get_cells_ephys_df(manifest_file='/local1/storage/celltypes_data/manifest.json')

for cell_id in cells_list[:]:
    # dynamics_path = os.path.join(params_dir, '{}.json'.format(cell_id))
    # node_dict = {'dynamics_params': dynamics_path,
    #                 'cell_name': "{}_stage1".format(cell_id)}
    dynamics_path = os.path.join(params_dir, 'optim_param_{}.json'.format(cell_id))
    dp = BiophysParams.from_json(dynamics_path)
    ih_params = {prop: ih_factors*float(dp.get_nested(prop)) for prop in props}
    ih_params['ih_factor'] = ih_factors
    path = sim_path.format(cell=cell_id, sim='psp')
    sm = sc.build_epsp_batch(cell_id, path, overwrite=overwrite,
                          inh=inh, n_duplicates=n_duplicates, transient_time=transient_time,
                          edge_dict=edge_dict, node_dict={}, linked_dict=ih_params)
    
    # account for junction potential
    vrest = cells_df.loc[cell_id, 'vrest'] - 14
    config_dict = {"conditions":{"v_init":vrest}}
    sm.config.update_nested(config_dict)
    sm.config.save()