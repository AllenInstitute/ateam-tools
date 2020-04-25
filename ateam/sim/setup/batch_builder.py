"""Module for building simulations consisting of multiple copies of a single cell,
with parameter variations, and keeping track of the cell parameters and outputs."""

from bmtk.builder.networks import NetworkBuilder
import collections
# from bmtk.simulator import bionet
# import bmtk.utils.sim_setup as setup
# import bmtk.utils.spike_trains.spikes_csv as spikes
import numpy as np
import pandas as pd
import h5py as h5
import os.path
from six import string_types
from collections import Iterable

def get_node_ids(net):
    node_ids = [node.node_id for node in net.nodes_iter()]
    return node_ids

def lookup_by_target(src, trg, prop_dict):
    return prop_dict.get(trg.node_id)

def build_props_combinatorial(props_base=None, n_duplicates=1, linked_dicts=None, **props_indep):
    """Build a props dict containing all permutations of a number of independent props,
    supplied as keyword args that each specifies a list of values (pass dict using **prop_dict)
    All combinations of independent props are formed, while properties specified together in linked_dicts
    are restricted to vary together (eg. distance_range min and max)
    If props_base dict is supplied, output will be a merged dict of base props and varied props.
    """
    props = props_base.copy() if props_base else {}
    # TODO: remove these lines, check for duplicate props
    keys_scalar = (key for key in props_indep.keys() if np.isscalar(props_indep[key]))
    props.update((key, props_indep.pop(key)) for key in keys_scalar)

    props_indep = collections.OrderedDict(props_indep)
    props_index_shape = [len(vals) for vals in props_indep.values()]

    # Validate linked props and add to the end of the index 
    linked_dicts = linked_dicts or []
    for props_linked in linked_dicts:
        n_linked = [len(prop) for prop in props_linked.values()]
        assert all(n_linked[0] == ni for ni in n_linked)
        n_linked = n_linked[0]
        props_index_shape += [n_linked]

    indices = np.indices(props_index_shape)    

    # Add independent props
    for i, prop in enumerate(props_indep.keys()):
        prop_array = np.array(props_indep[prop])
        props[prop] = np.repeat(prop_array[indices[i].flat], n_duplicates)

    # Add linked props - after independent, so can overwrite (not ideal, creates duplicates)
    n_props = len(props_indep)
    for j, props_linked in enumerate(linked_dicts):
        i = n_props + j
        for prop, vals in props_linked.items():
             prop_array = np.array(vals)
             props[prop] = np.repeat(prop_array[indices[i].flat], n_duplicates)

    props['N'] = int(n_duplicates*np.prod(props_index_shape))
    return props

def build_batch_node_props(node_props_base, n_duplicates=1, net_name='batch', linked_dicts=None, **vary_props):
    net = NetworkBuilder(net_name)
    node_props = build_props_combinatorial(node_props_base, n_duplicates, linked_dicts=linked_dicts, **vary_props)
    net.add_nodes(**node_props)
    net.build()
    return net

def build_input_net_simple(N=1, name='input', **props):
    net = NetworkBuilder(name)
    if 'num_input' in props.keys():
        num_input = props.pop('num_input')
        for k in range(num_input):
            net.add_nodes(N, model_type='virtual', pop_name = 'input_%s'%k, **props)
            
    else:
         net.add_nodes(N, model_type='virtual', **props)
    net.build()
    return net

def build_batch_edge_props(input_net, node_props_base, edge_props_base, n_duplicates=1, net_name='batch', linked_dicts=None, **vary_edge_props):
    net = NetworkBuilder(net_name)
    edge_props = build_props_combinatorial(n_duplicates=n_duplicates, linked_dicts=linked_dicts, **vary_edge_props)
    
    # Attach all props to the nodes just for record-keeping
    node_props = node_props_base.copy()
    node_props.update(edge_props)
    net.add_nodes(**node_props)
    
    cm = net.add_edges(source=input_net.nodes(), target=net.nodes(), **edge_props_base)

    N = edge_props.pop('N')
    node_ids = get_node_ids(net)
    for key, values in edge_props.items():
        prop_dict = dict(zip(node_ids, values))
        cm.add_properties(key, rule=lookup_by_target, rule_params={'prop_dict': prop_dict}, dtypes=values.dtype)
    net.build()
    return net

def split_dict_lists(full_dict):
    """Split a dict into separate dicts for singleton and list-like elements.
    Strings are considered as singletons, while any other iterable is considered a list.
    """
    
    list_dict = {}; single_dict = {}
    for key, val in full_dict.items():
        if isinstance(val, string_types) or not isinstance(val, Iterable):
            single_dict[key] = val
        else:
            list_dict[key] = val
    return single_dict, list_dict


def build_batch_all(sm, node_props, edge_props, input_props, n_duplicates=1,
                net_name='batch', linked_dicts=None, use_abs_paths=False,
                coordinate_props=None):
    abs_path_input = True if use_abs_paths==True or use_abs_paths=='input' else False
    abs_path_network = True if use_abs_paths==True or use_abs_paths=='network' else False
    node_props_base, node_props_vary = split_dict_lists(node_props)
    edge_props_base, edge_props_vary = split_dict_lists(edge_props)
    input_props_base, input_props_vary = split_dict_lists(input_props)

    # TODO: check and distinguish duplicate keys between dicts
    # TODO: specify linked as list, not as separate dict, to keep identity
    node_props_vary.update(edge_props_vary)
    node_props_vary.update(input_props_vary)

    all_props = build_props_combinatorial(n_duplicates=n_duplicates, linked_dicts=linked_dicts, **node_props_vary)
    N = all_props.pop('N')

    # Attach all props to the nodes just for record-keeping
    net = NetworkBuilder(net_name)
    input_net = None
    # node_props.update(edge_props_base) fails from dynamics_params
    node_props.update(input_props_base)
    node_props.update(all_props)

    if coordinate_props:

        temp_props = node_props.copy()
        for temp_props_key,temp_props_val in temp_props.items():
            if temp_props_key in node_props_vary:
                temp_props.pop(temp_props_key,None)
        for jj in range(N):

            temp_props.update({'x' : [coordinate_props['x'][jj]],
                               'y' : [coordinate_props['y'][jj]],
                               'z' : [coordinate_props['z'][jj]]})
            net.add_nodes(**temp_props)

    else:
        net.add_nodes(N=N, **node_props)

    sm.add_network(net)
    node_ids = get_node_ids(net)

    # For input props, combine base and varying
    input_props.update( (key, all_props[key]) for key in input_props_vary.keys() )
    rates = input_props.pop('input_rate', None)
    spike_time = input_props.pop('spike_time', None)
    spike_inputs = rates is not None or spike_time is not None
    current_inputs = 'amp' in input_props

    if spike_inputs:
        input_net = build_input_net_simple(N=N, **input_props)
        sm.add_network(input_net)

        multiple_inputs = False
        if rates is not None:
            if 'num_input' in input_props:
                multiple_inputs = True
                num_input = input_props.pop('num_input')
                rates = np.kron(np.ones(num_input),rates)
            sm.add_spike_input_poisson(input_net.name, rates, use_abs_paths=abs_path_input, **input_props)

        if spike_time is not None:
            sm.add_spike_input_vector(input_net.name, [spike_time], use_abs_paths=abs_path_input)

        # For edge props, keep base and varying separate
        edge_props_vary.update( (key, all_props[key]) for key in edge_props_vary.keys() )
        if multiple_inputs:
            connections = []
            for ii in range(num_input):
                connections.append( net.add_edges(source=input_net.nodes(pop_name='input_%s'%ii),
                                target=net.nodes(), iterator = 'paired',  **edge_props_base))
        else:
            connections = [net.add_edges(source=input_net.nodes(), target=net.nodes(), iterator='paired', **edge_props_base)]

        for cm in connections:
            for key, values in edge_props_vary.items():
                prop_dict = dict(zip(node_ids, values))
                cm.add_properties(key, rule=lookup_by_target, rule_params={'prop_dict': prop_dict}, dtypes=values.dtype)
    
    if current_inputs:
        sm.add_current_clamp_complex(net.name, input_props, use_abs_paths=abs_path_input)

    net.build()
    sm.save_network_files(use_abs_paths=abs_path_network)
    return net, input_net

def read_node_props_batch(net_folder_path):
    # simple single-pop network, no need for fancy
    net = "batch"
    nodes_h5_file = "{net}_nodes.h5".format(net=net)
    nodes_h5_path = os.path.join(net_folder_path, nodes_h5_file)
    return read_node_props(nodes_h5_path, net)

def read_node_props(nodes_h5_path, net, nodeset=0):
    # Could also use analyzer.nodes_table or merged node_table, or sonata api?
    group_path = "nodes/{net}/{nodeset}".format(net=net, nodeset=nodeset)
    with h5.File(nodes_h5_path) as nodes_h5:
        node_group = nodes_h5[group_path]
        nodes_df = pd.DataFrame.from_dict(dict(node_group.iteritems()))
    return nodes_df