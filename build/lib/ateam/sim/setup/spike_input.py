import pandas as pd
import numpy as np
import csv
import h5py
from six import string_types
from bmtk.builder.node_pool import NodePool

class SpikeInput(object):
    def __init__(self, nodes, populations=None):
        """SpikeInput class creates and stores input spike trains to a set of nodes.
        nodes can be either the path to a nodes h5 file, or a NodePool object.
        """
        try:
            if isinstance(nodes, string_types):
                nodes_h5 = h5py.File(nodes, 'r')
                nodes_grp = nodes_h5['/nodes']
                if populations is None:
                    populations = nodes_grp.keys()

                # TODO: Need a way to Use sonata library without having to use node-types
                node_ids = []
                for node_pop in populations:
                    node_ids.extend(nodes_grp[node_pop]['node_id'])
                self._network = nodes.split('_')[0]

            else:
            # isinstance(nodes, NodePool):
                node_ids = [node.node_id for node in nodes]
                self._network = nodes.network_name
        except:
            Exception("Can't read nodes.")
        self._spikemap = {n: [] for n in nodes}
        self._nodes = node_ids

    @property
    def network(self):
        return self._network

    def set_times_all(self, spike_times):
        self._spikemap = {n: spike_times for n in self._nodes}

    def save_csv(self, csv_file_name, in_ms=True):
        conv = 1000.0 if in_ms else 1.0

        with open(csv_file_name, 'w') as csv_file:
            csv_writer = csv.writer(csv_file, delimiter=' ')
            csv_writer.writerow(['gid', 'spike-times'])
            for gid, rate_gen in self._spikemap.items():
                csv_writer.writerow([gid, ','.join(str(r*conv) for r in rate_gen)])