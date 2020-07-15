import ateam.sim.singlecell as sc
import os
import os.path

cells_path = "/allen/aibs/mat/ateam_shared/Human_Model_Fit_Metrics"
cells_all = [int(name) for name in os.listdir(cells_path) if name.isdigit()]

# bootstrap from previous run
cells_path = "/allen/aibs/mat/tmchartrand/bmtk_networks/singlecell/round1"
cells_inh = os.listdir(os.path.join(cells_path,"IN"))
cells_exc = os.listdir(os.path.join(cells_path,"PC"))

base_path = "/allen/aibs/mat/tmchartrand/bmtk_networks/singlecell/round3/IN"
sc.build_singlecell_sims(cells_inh, base_path, inh=True, overwrite=True, level='network')

base_path = "/allen/aibs/mat/tmchartrand/bmtk_networks/singlecell/round3/PC"
sc.build_singlecell_sims(cells_exc, base_path, overwrite=True, level='network')