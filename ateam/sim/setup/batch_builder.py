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

def get_node_ids(net):
    node_ids = [node.node_id for node in net.nodes_iter()]
    return node_ids

def lookup_by_target(src, trg, prop_dict):
    return prop_dict.get(trg.node_id)

def build_props_combinatorial(props_base=None, n_duplicates=1, props_linked=None, **props_indep):
    """Build a props dict containing all permutations of a number of independent props
    Each property takes values from its own set.
    """
    props = props_base.copy() if props_base else {}
    props_dict = collections.OrderedDict(props_indep)
    props_shape = [len(vals) for vals in props_dict.values()]
    indices = np.indices(props_shape)

    if props_linked:
        n_linked = [len(prop) for prop in props_linked.values()]
        assert all(n_linked[0] == ni for ni in n_linked)
        n_linked = n_linked[0]
        # use the final dimension (-1) for linked props
        props_shape = props_shape + [n_linked]
        indices = np.indices(props_shape)
        for prop, vals in props_linked.items():
             prop_array = np.array(vals)
             props[prop] = np.repeat(prop_array[indices[-1].flat], n_duplicates)

    for i, prop in enumerate(props_dict.keys()):
        prop_array = np.array(props_dict[prop])
        props[prop] = np.repeat(prop_array[indices[i].flat], n_duplicates)

    props['N'] = n_duplicates*np.prod(props_shape)
    return props

def build_batch_node_props(node_props_base, n_duplicates=1, net_name='batch', props_linked=None, **vary_props):
    net = NetworkBuilder(net_name)
    node_props = build_props_combinatorial(node_props_base, n_duplicates, props_linked=props_linked, **vary_props)
    net.add_nodes(**node_props)
    net.build()
    return net

def build_input_net_simple(N=1):
    net = NetworkBuilder('input')
    net.add_nodes(N, model_type='virtual')
    net.build()
    return net

def build_batch_edge_props(input_net, node_props_base, edge_props_base, n_duplicates=1, net_name='batch', props_linked=None, **vary_edge_props):
    net = NetworkBuilder(net_name)
    edge_props = build_props_combinatorial(n_duplicates=n_duplicates, props_linked=props_linked, **vary_edge_props)
    # Attach all props to the nodes just for record-keeping
    node_props = node_props_base.copy()
    node_props.update(edge_props)
    net.add_nodes(**node_props)
    
    cm = net.add_edges(source=input_net.nodes(), target=net.nodes(), **edge_props_base)

    N = edge_props.pop('N')
    node_ids = get_node_ids(net)
    for key, values in edge_props.items():
        prop_dict = dict(zip(node_ids, values))
        cm.add_properties(key, rule=lookup_by_target, rule_params={'prop_dict': prop_dict})
    net.build()
    return net

def read_node_props_batch(net_folder_path):
    net = "batch"
    nodes_h5_file = "{net}_nodes.h5".format(net=net)
    nodes_h5_path = os.path.join(net_folder_path, nodes_h5_file)
    return read_node_props(nodes_h5_path, net)

def read_node_props(nodes_h5_path, net, nodeset=0):
    group_path = "nodes/{net}/{nodeset}".format(net=net, nodeset=nodeset)
    with h5.File(nodes_h5_path) as nodes_h5:
        node_group = nodes_h5[group_path]
        nodes_df = pd.DataFrame.from_dict(dict(node_group.iteritems()))
    return nodes_df