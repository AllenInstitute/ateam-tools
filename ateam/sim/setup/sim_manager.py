import os
import shutil
import subprocess
import itertools
import csv
from collections import defaultdict
# import bmtk.simulator.utils.config as config
from bmtk.simulator.utils.config import ConfigDict
from bmtk.utils.io.spike_trains import PoissonSpikesGenerator
import bmtk.builder.networks as buildnet
import bmtk.utils.sim_setup as setup
import bmtk.analyzer.visualization.spikes as vs

from .config_class import ConfigBuilder
from .spike_input import SpikeInput
from ateam.sim.run import runner

ConfigClass = ConfigBuilder

class SimManager(object):
    def __init__(self, config_path="./config.json", sim_folder='', import_nodes=False):
        # TODO: relax constraint on sim_folder as parent dir?
        """Create a SimManager for a simulation defined by a config_file,
        using the parent folder as the simulation folder.
        If config file is not specified, looks for config.json in current directory.
        """
        config_path = os.path.abspath(os.path.join(sim_folder, config_path))
        # assert os.path.isfile(config_path)
        self.config = ConfigClass.from_json(config_path)
        self.configdict = ConfigDict.load(config_path)
        self.sim_folder = os.path.dirname(config_path)
        self._networks_active = {}
        # Dict of net_name: nodes file dict
        self._nodes_dict = {}
        # Dict of (src_name, trg_name): list of edge file dicts
        self._edges_dict = defaultdict(list)
        self.load_networks(import_nodes=import_nodes)

    @classmethod
    def from_template(cls, config_template, config_file="config.json", sim_folder=None, config_path=None, overwrite=False):
        """
        Create a SimManager from template config file in a new simulation folder.
        Creates folder if it doesn't exist.
        """
        if config_path:
            config_path = os.path.expandvars(os.path.expanduser(config_path))
            sim_folder = os.path.dirname(config_path)
        else:
            sim_folder = sim_folder or os.getcwd()
            sim_folder = os.path.expandvars(os.path.expanduser(sim_folder))
            config_path = os.path.join(sim_folder, config_file)
        # TODO: work with relative paths?
        if not os.path.exists(sim_folder):
            os.makedirs(sim_folder)

        if overwrite or not os.path.isfile(config_path):
            ConfigClass.load_template(config_template, config_path, shared_components=True)
        else:
            Warning("Config file already exists: loading config without template.")

        return cls(config_path)

    
    @property
    def config_path(self):
        return self.config.path

    @property
    def network_builders(self):
        return self._networks_active

    @property
    def networks(self):
        all_nets = list(self._networks_active.keys()) + list(self._nodes_dict.keys())
        return all_nets

    @property
    def networks_saved(self):
        return self._nodes_dict.keys()

    @property
    def files_list(self):
        # return itertools.chain(*self.files_dict.values())
        raise NotImplementedError()

    def new_network(self, net_name):
        net = buildnet.NetworkBuilder(net_name)
        self.add_network(net)
        return net

    def add_network(self, network):
        # TODO: check for name collision
        self._networks_active[network.name] = network

    def add_networks(self, networks):
        for network in networks:
            self.add_network(network)

    def clear_output(self):
        out_dir = os.path.join(self.sim_folder, "output")
        shutil.rmtree(out_dir)

    def save_network_files(self, net_folder_name='', use_abs_paths=False):
        net_path = self.abspath(net_folder_name)
        if use_abs_paths:
            trim_path = lambda path: path
        else:
            trim_path = lambda path: os.path.relpath(path, self.sim_folder)
            
        for name, net in self._networks_active.items():
            if name not in self.networks_saved: # doesn't try to save already loaded networks
                nodes_dict, edge_dicts = net.save(net_path)
                # Convert back to relative paths for config
                for key, path in nodes_dict.items():
                    nodes_dict[key] = trim_path(path)
                self._nodes_dict.update({name: nodes_dict})
                for edgeset in edge_dicts:
                    for key, path in edgeset.items():
                        edgeset[key] = trim_path(path)
                    self._edges_dict[edges_net_pair(edgeset)].append(edgeset)

        nodes = self._nodes_dict.values()
        edges = sum(self._edges_dict.values(), [])
        self.config.update(networks={'nodes':nodes, 'edges':edges})
        self.config.save()

    def save_complete(self, folder_path=''):
        """Save self-contained folder with all simulation files."""
        # TODO: could just copy folder as a first step, but really need to identify external files
        folder_path = folder_path or self.sim_folder
        raise NotImplementedError()

    def save_copy(self, folder_path, overwrite=False):
        """Save a copy of the current sim folder contents to a new directory.
        Remove the output directory from the new folder."""
        if os.path.exists(folder_path):
            if overwrite:
                shutil.rmtree(folder_path)
            else:
                Warning("Folder exists; aborting.")
                return
        shutil.copytree(self.sim_folder, folder_path)
        config_name = os.path.basename(self.config_path)
        sm = SimManager(config_name, sim_folder=folder_path)
        sm.clear_output()
        return sm

    def load_networks(self, import_nodes=False):
        """Loads file paths for networks specified in the config"""
        # Need to use configdict to resolve paths
        if self.configdict.with_networks:
            nets_dict = self.configdict.networks
            # gets network block from config file
            nodes_raw = nets_dict['nodes']
            edges_raw = nets_dict['edges']
            if import_nodes:
                nodes = {}
                for nodeset in nodes_raw:
                    name = nodes_net_name(nodeset)
                    net = self.new_network(name)
                    net.import_nodes(nodeset['nodes_file'], nodeset['node_types_file'], population=None)
                    
            nodes = {nodes_net_name(nodeset):nodeset for nodeset in nodes_raw}
            self._nodes_dict.update(nodes)

            for edgeset in edges_raw:
                self._edges_dict[edges_net_pair(edgeset)].append(edgeset)
        
    def update_node_type_props(self, net_name, props):
        assert(net_name in self.networks_saved)
        node_types_file = self._nodes_dict[net_name]['node_types_file']
        update_csv(self.abspath(node_types_file), props)

    def update_edge_type_props(self, src_name, dest_name, props):
        edges = self._edges_dict.get((src_name, dest_name))
        if not edges:
            raise Exception("Could not find edges for the specified networks.")
        if len(edges) > 1:
            Warning("Multiple edge files exist for given pair of networks.")
        edge_dict = edges[0]
        update_csv(self.abspath(edge_dict['edge_types_file']), props)


### Configure inputs and modules
    def abspath(self, filename):
        return os.path.join(self.sim_folder, filename)

    def add_spike_input(self, input_file, net_name, trial=None):
        """Add specified spikeinput file to the config.
        Note that node ids in the file must match those for the specified net.
        """
        # assert(self._networks_active.has_key(net_name))
        ext = os.path.splitext(input_file)[1][1:]
        inputs = {net_name: 
            {
            'input_type': 'spikes',
            'module': ext,
            'input_file': input_file,
            'node_set': net_name,
            'trial': trial
            }}
        self.config.update_nested(inputs=inputs)
    
    def write_spikeinput_vector(self, net_name, times, spike_file_name='spike_input.csv',use_abs_paths=False):
        """Write a new spikeinput file from a vector of times and add it to the config.
        All cells are assigned the same spike times."""
        spikes = SpikeInput(self._networks_active[net_name].nodes())
        spikes.set_times_all(times)
        if use_abs_paths:
            spike_file = os.path(spike_file_name)
        else:    
            spike_file = self.abspath(spike_file_name)
        spikes.save_csv(spike_file)
        self.add_spike_input(spike_file_name, net_name)

    def write_spikeinput_poisson(self, net_name, rate, tstart = 0, tstop=2000,\
                                 spike_file_name='spike_input.h5',use_abs_paths=False):
        """Write a new spikeinput file for independent Poisson spiking and add it to the config."""
        net = self._networks_active[net_name]
        node_ids = [node.node_id for node in net.nodes_iter()]
        psg = PoissonSpikesGenerator(node_ids, rate, tstart = tstart, tstop=tstop)
        spike_file = self.abspath(spike_file_name)
        
        if use_abs_paths:
            spike_file = os.path.abspath(self.abspath(spike_file_name))
            spike_file_name = os.path.abspath(self.abspath(spike_file_name))
            
        psg.to_hdf5(spike_file)
        self.add_spike_input(spike_file_name, net_name)

    def add_ecp_report(self, electrode_file=None, cells='all', file_name='ecp.h5', locs=[[0,0,0]], save_contributions=False):
        if electrode_file is None:
            electrode_file = 'electrode.csv'
            self.write_electrode_file(locs, self.abspath(electrode_file))

        reports = {'ecp_report': {
            'cells': cells,
            'variable_name': 'v',
            'module': 'extracellular',
            'electrode_positions': electrode_file,
            'file_name': file_name,
            'electrode_channels': 'all',
            'contributions_dir': 'ecp_contributions' if save_contributions else None
        }}
        self.config.update_nested(reports=reports)


    @staticmethod
    def write_electrode_file(locs, csv_file_name):
        import csv
        with open(csv_file_name, 'w') as csv_file:
            csv_writer = csv.writer(csv_file, delimiter=' ')
            csv_writer.writerow(['channel', 'x_pos', 'y_pos', 'z_pos'])
            for i, loc in enumerate(locs):
                csv_writer.writerow([i] + [str(x) for x in loc])


    def add_membrane_report(self, name='membrane_report', variables=['v'], cells='all', sections='soma',
                            file_name=None, **kwargs):
        reports = {name: {
            'module': 'membrane_report',
            'variable_name': variables,
            'cells': cells,
            'file_name': file_name or '{name}.h5'.format(name=name),
            'sections': sections
            }}
        reports[name].update(**kwargs)
        self.config.update_nested(reports=reports)
        self.config.save()
        
    def nodes_file(self, net):
        # TODO: check net in networks first?
        return self.abspath("{}_nodes.h5".format(net))
        
    def node_types_file(self, net):
        # TODO: check net in networks first?
        return self.abspath("{}_node_types.csv".format(net))

    @property
    def spikes_file(self):
        return self.configdict.spikes_file

    @property
    def sim_time(self):
        return self.config['run']['tstop']

    @sim_time.setter
    def sim_time(self, time):
        self.config.update_nested(run={"tstop": time})

    # TODO: remove below
    def set_sim_time(self, time):
        self.sim_time = time

    def run_bionet(self):
        self.config.save()
        return runner.run_bionet(self.config_path)
        
    def run_bionet_mpi(self, ncores=1):
        self.config.save()
        runner.run_bionet_mpi(self.config_path, ncores)

    def plot_raster(self, net, **kwargs):
        return vs.plot_spikes(self.nodes_file(net), self.node_types_file(net), self.spikes_file, **kwargs)

    def plot_rates(self, net, **kwargs):
        return vs.plot_rates(self.nodes_file(net), self.node_types_file(net), self.spikes_file, **kwargs)


def nodes_net_name(nodeset):
    """Extract the network name from nodes file dict as specified in config"""
    nodes_file = nodeset['nodes_file']
    net_name = os.path.basename(nodes_file).split('_')[0]
    return net_name

def edges_net_pair(edgeset):
    """Extract the source and target network names from edges file dict as specified in config.
    Return a tuple (source name, target name)."""
    edge_file = edgeset['edges_file']
    src, dest = os.path.basename(edge_file).split('_')[0:2]
    return (src, dest)

def update_csv(csv_path, props):
    with open(csv_path) as csvfile:
        reader = csv.DictReader(csvfile, delimiter=' ')
        rows = [dict(row, **props) for row in reader]
    with open(csv_path, 'w',) as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=rows[0].keys(), delimiter=' ')
        writer.writeheader()
        writer.writerows(rows)