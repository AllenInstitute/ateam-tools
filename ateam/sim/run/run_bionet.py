# -*- coding: utf-8 -*-

"""Simulates an example network of 14 cell receiving two kinds of exernal input as defined in configuration file"""

import os, sys
from bmtk.simulator import bionet
# Register our modified axon processing with bionet (call after import bionet)
from ateam.sim import cell_functions 

def create_sim(config_file):
    conf = bionet.Config.from_json(config_file, validate=True)
    conf.build_env()

    graph = bionet.BioNetwork.from_config(conf)
    sim = bionet.BioSimulator.from_config(conf, network=graph)
    return sim

def run(config_file):
    conf = bionet.Config.from_json(config_file, validate=True)
    conf.build_env()

    graph = bionet.BioNetwork.from_config(conf)
    sim = bionet.BioSimulator.from_config(conf, network=graph)
    sim.run()
    bionet.nrn.quit_execution()

def main():
    if __file__ != sys.argv[-1]:
        run(sys.argv[-1])
    else:
        run('config.json')

if __name__ == '__main__':
    main()
