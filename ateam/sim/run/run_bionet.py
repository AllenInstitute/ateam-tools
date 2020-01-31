# -*- coding: utf-8 -*-

import os, sys
from bmtk.simulator import bionet
# Register our modified axon processing with bionet (call after import bionet)
from ateam.sim import cell_functions 

def create_sim(config_file, setup_output=True):
    conf = bionet.Config.from_json(config_file, validate=True)
    conf.build_env(setup_output=setup_output)

    graph = bionet.BioNetwork.from_config(conf)
    sim = bionet.BioSimulator.from_config(conf, network=graph)
    return sim

def run(config_file):
    conf = bionet.Config.from_json(config_file, validate=True)
    conf.build_env()

    graph = bionet.BioNetwork.from_config(conf)
    sim = bionet.BioSimulator.from_config(conf, network=graph)
    sim.run()
    # bionet.nrn.quit_execution()
    bionet.nrn.reset()

def main():
    if __file__ != sys.argv[-1]:
        run(sys.argv[-1])
    else:
        run('config.json')

if __name__ == '__main__':
    main()
