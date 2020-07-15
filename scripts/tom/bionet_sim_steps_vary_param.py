import ateam.sim.setup.default_props as defaults
from ateam.sim.setup import SimManager
from ateam.sim.run import run_bionet
import numpy as np
from ateam.analysis.bmtk.cell_vars import plot_v
import os
import os.path
from bmtk.simulator.bionet.biophys_params import BiophysParams

# Specify the sim folder and a template file
temp_dir = '/local1/storage/temp/bmtk/'
os.chdir(temp_dir)
new_config = temp_dir+'batch_config.json'
config_template = "/allen/aibs/mat/tmchartrand/bmtk_networks/biophys_components_shared/default_config.json"

sm = SimManager.from_template(config_template=config_template, overwrite=True, config_path=new_config)

cell_id = 525133308
# node_props = defaults.cellprops_active(cell_id)
params_path = "absolute/path/to/params.json"
node_props = {
        'cell_name': cell_name(cell_id),
        'morphology': morph_file(cell_id),
        'dynamics_params': params_path,
        'model_type': 'biophysical',
        'model_template': 'ctdb:Biophys1.hoc',
        'model_processing': 'aibs_allactive_ani_directed',
    }
dynamics_params = BiophysParams.from_json(params_path)
net = sm.new_network('single')

N = 3
prop = 'genome.soma.cm'
scale_factor = np.linspace(0.5, 2, N)
node_props[prop] = scale_factor*float(dynamics_params.get_nested(prop))
net.add_nodes(N=N, **node_props)

sm.sim_time = 100
sm.config.update_nested(inputs={
     "current_clamp": {
      "input_type": "current_clamp",
      "module": "IClamp",
      "node_set": "all",
      "amp": 0.120,
      "delay": 5.0,
      "duration": 1000.0
    }})

sm.add_membrane_report()
sm.save_network_files()
run_bionet.run(sm.config_path)
plot_v(sm.config_path, show=True)