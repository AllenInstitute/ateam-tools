import ateam.sim.singlecell_analysis as sca
import ateam.analysis.dataframes as dfa
import pandas as pd
import os
import os.path
from functools import partial

def fit_path(base_path, cells_list=None):
    cells_list = cells_list or os.listdir(base_path)
#     Need to include morphology column to track cell id for hof
    fit_function = partial(sca.fit_df_all, base_path=base_path)
    fit_all = map(fit_function, cells_list)
    cells_fit_df = pd.concat(fit_all, axis=0, keys=cells_list, names=['cell'])
    
    cells_fit_df = dfa.flatten_columns(cells_fit_df).reset_index(level="target_sections")
    return cells_fit_df

path_all = "/allen/aibs/mat/tmchartrand/bmtk_networks/singlecell/round3"

cells_fit_df = pd.concat([sca.fit_path(os.path.join(path_all, 'IN')), sca.fit_path(os.path.join(path_all, 'PC'))])
cells_fit_df['specimen_id'] = [int(cell) for cell in cells_fit_df.index.values]
