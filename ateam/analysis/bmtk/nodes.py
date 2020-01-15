import pandas as pd
import matplotlib.pyplot as plt
import h5py
from ateam.sim.setup import SimManager

def create_node_table(node_file, node_type_file, group_key=None, exclude=[], node_props=None):
    """Creates a merged nodes and node_types dataframe with excluded items removed. Returns a dataframe.
    For properties in the nodes file, includes all by default, or a subset specified in node_props (as a list of names).
    Forked from bmtk.analyzer.visualization.spikes, """
    node_types_df = pd.read_csv(node_type_file, sep=' ', index_col='node_type_id')
    # TODO: use bmtk classes here!
    with h5py.File(node_file, 'r') as nodes_h5:
        node_pop_name = nodes_h5['/nodes'].keys()[0]
        nodes_grp = nodes_h5['/nodes'][node_pop_name]
        nodes_df = pd.DataFrame({key: nodes_grp[key] for key in ['node_id', 'node_type_id']})
        if node_props is None:
            keys = nodes_grp['0'].keys()
            props_dict = {key: nodes_grp['0'][key] for key in keys if key != 'dynamics_params'}
            nodes_df = nodes_df.join(pd.DataFrame.from_dict(props_dict))
        else:
            for key in node_props:
                nodes_df[key] = nodes_grp['0'][key]
        nodes_df.set_index('node_id', inplace=True)

    full_df = pd.merge(left=nodes_df, right=node_types_df, how='left', left_on='node_type_id', right_index=True)

    if group_key is not None and len(exclude) > 0:
        # Make sure sure we group-key exists as column
        if group_key not in full_df:
            raise Exception('Could not find column {}'.format(group_key))
        full_df = full_df[~full_df[group_key].isin(exclude)]

    return full_df

def plot_soma_positions(config_file, netname, group_key, color_dict):
    sm = SimManager(config_file)
    nodes_df = create_node_table(sm.nodes_file(netname), sm.node_types_file(netname), group_key=group_key, node_props=['positions'])
    groups = nodes_df.groupby(group_key)

    for group_name, group_df in groups:
        positions = group_df['positions']
        plt.scatter([pos[0] for pos in positions], [pos[1] for pos in positions], color=color_dict.get(group_name, 'k'))
    
    plt.xlabel("X pos ($\mu$m)")
    plt.ylabel("Y pos ($\mu$m)")