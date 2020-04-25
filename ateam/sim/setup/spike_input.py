import pandas as pd
import numpy as np
import csv
import h5py
from six import string_types
from bmtk.builder.node_pool import NodePool

class NodeInput(object):
    def __init__(self, nodes, populations=None):
        """Creates and stores inputs to a set of nodes.
        nodes can be either the path to a nodes h5 file, or a NodePool object.
        Based on bmtk.utils.spike_trains.spikes_csv.SpikesGenerator
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
                self.network = nodes.split('_')[0]

            else:
            # isinstance(nodes, NodePool):
                node_ids = [node.node_id for node in nodes]
                self.network = nodes.network_name
        except:
            raise Exception("Can't read nodes.")
        self.nodes = node_ids
        self._spike_inputs = None
        self._current_inputs = None

    def set_current_inputs(self, input_dict, csv_file_name):
        self._current_inputs = input_dict
        input_dict['gid'] = self.nodes
        current_data = pd.DataFrame.from_dict(input_dict)
        current_data.to_csv(csv_file_name, index=False)


    def set_spike_inputs_all_nodes(self, spike_times, csv_file_name):
        """Assign the same list of spike times to all nodes.
        Spike times in ms.
        """
        self._spike_inputs = {n: spike_times for n in self.nodes}
        self.save_spikes_csv(csv_file_name)

    def save_spikes_csv(self, csv_file_name):
        with open(csv_file_name, 'w') as csv_file:
            csv_writer = csv.writer(csv_file, delimiter=' ')
            csv_writer.writerow(['gid', 'spike-times'])
            for gid, rate_gen in self._spike_inputs.items():
                csv_writer.writerow([gid, ','.join(str(r*conv) for r in rate_gen)])