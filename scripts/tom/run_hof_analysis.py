import ateam.sim.singlecell_analysis as sca
import ateam.analysis.dataframes as dfa

import pandas as pd
import os
import multiprocessing as mp
from functools import partial
from itertools import chain

base_path = "/allen/aibs/mat/tmchartrand/bmtk_networks/singlecell/round3/hof/"


# cells_dirs = [os.path.join(hof_path, f) for f in os.listdir(hof_path)]
# models_list = chain([os.path.join(cell,model) for model in os.listdir(cell)] for cell in cells_dirs)

hof_range = range(10)
cells_list = os.listdir(base_path)
def fit_cell(cell): 
    return [sca.fit_df_all(cell, base_path=base_path, extra=["cell_id", "cell_name"], hof_id=hof_id) for hof_id in hof_range]

p = mp.Pool()
print("{} processors".format(mp.cpu_count()))
print("{} workers".format(p._processes))
fit_all = p.imap(fit_cell, cells_list)
# fit_all = map(fit_cell, cells_list)

cells_fit_df = pd.concat(chain(*fit_all), axis=0, sort=True)
cells_fit_df = dfa.flatten_columns(cells_fit_df).reset_index(level="target_sections")

cells_fit_df.to_csv(base_path + "../hof_round3_pc_2.csv")