#!/usr/bin/env python

import argparse
import os
import ateam.sim.singlecell as singlecell

parser = argparse.ArgumentParser()
parser.add_argument("base_path", default=".")
args = parser.parse_args()
cells_list = os.listdir(args.base_path)
singlecell.submit_singlecell_sims(cells_list, args.base_path, rerun=False)