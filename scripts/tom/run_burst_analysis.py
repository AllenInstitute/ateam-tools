import ateam.analysis.ephys_features.burst as burst
import ateam.analysis.dataframes as dfa

import pandas as pd
import multiprocessing as mp
from functools import partial
from itertools import chain

cells_csv = "pv_human_jim_tri.csv"
index_col="id"
out_path = "burst_analysis_04-03"

p = mp.Pool()
print("{} processors".format(mp.cpu_count))
print("{} workers".format(p._processes))
jim_df = pd.read_csv(cells_csv, index_col=index_col)
specimen_ids = jim_df.index.values
results = p.map(partial(burst.process_cell, passed_only=False, use_silence=True), specimen_ids)
sweep_results_df = pd.DataFrame.from_records(chain(*results))

sweep_results_df.to_csv(out_path + "_sweeps.csv")

sweep_results_df.drop(columns=["sweep_num"], inplace=True)
cell_results_df = dfa.flatten_columns(sweep_results_df.groupby("cell_id").agg(['mean','max']))
cell_results_df.to_csv(out_path + "_cells.csv")