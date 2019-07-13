#!/usr/bin/env python

import argparse
import os
import ateam.sim.singlecell as singlecell
import ateam.sim.run.run_hpc as run_hpc
import numpy as np

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


parser = argparse.ArgumentParser()
parser.add_argument("base_path", default=".")
args = parser.parse_args()
cells_list = os.listdir(args.base_path)
submit_singlecell_sims(cells_list, args.base_path, rerun=False)