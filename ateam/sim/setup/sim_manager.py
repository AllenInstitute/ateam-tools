import os
import subprocess
import itertools
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
    def __init__(self, config_path="./config.json"):
        """Create a SimManager for a simulation defined by a config_file,
        using the parent folder as the simulation folder.
        If config file is not specified, looks for config.json in current directory.
        """
        assert os.path.isfile(config_path)
        self.config = ConfigClass.from_json(config_path)
        self.sim_folder = os.path.dirname(config_path)
        # self.morphologies_dir = self.config.morphologies_dir
        # self.fit_dir = self.config.biophysical_neuron_models_dir
        self._networks = {}
        self._files_dict = {}

    @classmethod
    def from_template(cls, config_template, config_file="config.json", sim_folder=None, overwrite=False):
        """
        Create a SimManager from template config file in a new simulation folder.
        Creates folder if it doesn't exist.
        """

        sim_folder = sim_folder or os.getcwd()
        sim_folder = os.path.expandvars(os.path.expanduser(sim_folder))
        config_path = os.path.join(sim_folder, config_file)
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
    def networks(self):
        return self._networks.keys()

    @property
    def files_dict(self):
        return self._files_dict

    @property
    def files_list(self):
        return itertools.chain(*self.files_dict.values())


    def add_network(self, network):
        # TODO: check for name collision
        self._networks[network.name] = network

    def add_networks(self, networks):
        for network in networks:
            self.add_network(network)

    def save_network_files(self, net_folder_name=''):
        net_path = os.path.join(self.sim_folder, net_folder_name)


        for name, net in self._networks.items():
            if not self._files_dict.get(name):
                files = net.save(net_path)
                self._files_dict[name] = files
        # TODO: use Config class directly, not file
        # class could simplify using $NETWORK_DIR from manifest

        net_json = setup.get_network_block(files=self.files_list)
        self.config.update_nested(networks=net_json)
        self.config.save()


    def save_complete(self, folder_path=''):
        """Save self-contained folder with all simulation files."""
        folder_path = folder_path or self.sim_folder
        raise NotImplementedError()


    def load_networks(self):
        raise NotImplementedError()

### Configure inputs and modules
    def path(self, filename):
        return os.path.join(self.sim_folder, filename)

    def add_spike_input(self, input_file, net_name, trial=None):
        """Add specified spikeinput file to the config.
        Note that node ids in the file must match those for the specified net.
        """
        # assert(self._networks.has_key(net_name))
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
    
    def write_spikeinput_vector(self, net_name, times, spike_file_name='spike_input.csv'):
        """Write a new spikeinput file from a vector of times and add it to the config.
        All cells are assigned the same spike times."""
        spikes = SpikeInput(self._networks[net_name].nodes())
        spikes.set_times_all(times)
        spikes.save_csv(self.path(spike_file_name))
        self.add_spike_input(spike_file_name, net_name)

    def write_spikeinput_poisson(self, net_name, rate, tstop=2000, spike_file_name='spike_input.h5'):
        """Write a new spikeinput file for independent Poisson spiking and add it to the config."""
        net = self._networks[net_name]
        node_ids = [node.node_id for node in net.nodes_iter()]
        psg = PoissonSpikesGenerator(node_ids, rate, tstop=tstop)
        psg.to_hdf5(self.path(spike_file_name))
        self.add_spike_input(spike_file_name, net_name)

    def add_ecp_report(self, electrode_file=None, cells='all', file_name='ecp.h5'):
        if electrode_file is None:
            electrode_file = 'electrode.csv'
            self.write_electrode_file([[0,0,0]], electrode_file)

        reports = {'ecp_report': {
            'cells': cells,
            'variable_name': 'v',
            'module': 'extracellular',
            'electrode_positions': electrode_file,
            'file_name': file_name,
            'electrode_channels': 'all',
            'contributions_dir': 'ecp_contributions'
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


    def add_membrane_report(self, variables=['v'], cells='all', sections='soma', file_name='cell_vars.h5'):
        reports = {'membrane_report': {
            'module': 'membrane_report',
            'variable_name': variables,
            'cells': cells,
            'file_name': file_name,
            'sections': sections,
        }}
        self.config.update_nested(reports=reports)
        self.config.save()

    def nodes_file(self, net):
        # TODO: check net in networks first?
        return self.path("{}_nodes.h5".format(net))
        
    def node_types_file(self, net):
        # TODO: check net in networks first?
        return self.path("{}_node_types.csv".format(net))

    @property
    def spikes_file(self):
        # TODO: incorporate here not in ConfigDict?
        conf = ConfigDict.load(self.config.path)
        return conf.spikes_file

    def plot_raster(self, net, group_key=None):
        vs.plot_spikes(self.nodes_file(net), self.node_types_file(net), self.spikes_file, group_key=group_key)

    def plot_rates(self, net, group_key=None):
        vs.plot_rates(self.nodes_file(net), self.node_types_file(net), self.spikes_file, group_key=group_key)

    @property
    def sim_time(self):
        return self.config['run']['tstop']

    # TODO: remove
    def set_sim_time(self, time):
        self.sim_time = time
    
    @sim_time.setter
    def sim_time(self, time):
        self.config.update_nested(run={"tstop": time})

    def run_bionet(self):
        self.config.save()
        runner.run_bionet(self.config_path)
        
    def run_bionet_mpi(self, ncores=1):
        self.config.save()
        runner.run_bionet_mpi(self.config_path, ncores)
