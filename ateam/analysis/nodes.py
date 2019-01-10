import pandas as pd
import h5py

def create_node_table(node_file, node_type_file, group_key=None, exclude=[], node_props=None):
    """Creates a merged nodes and node_types dataframe with excluded items removed. Returns a dataframe.
    For properties in the nodes file, includes all by default, or a subset specified in node_props (as a list of names).
    Forked from bmtk.analyzer.visualization.spikes, """
    node_types_df = pd.read_csv(node_type_file, sep=' ', index_col='node_type_id')
    
    with h5py.File(node_file) as nodes_h5:
        node_pop_name = nodes_h5['/nodes'].keys()[0]
        nodes_grp = nodes_h5['/nodes'][node_pop_name]
        nodes_df = pd.DataFrame({key: nodes_grp[key] for key in ['node_id', 'node_type_id']})
        if node_props is None:
            nodes_df = nodes_df.join(pd.DataFrame.from_dict(dict(nodes_grp['0'].iteritems())))
        else:
            for key in node_props:
                nodes_df[key] = nodes_grp['0'][key]
        nodes_df.set_index('node_id', inplace=True)

    # nodes_df = pd.read_csv(node_file, sep=' ', index_col='node_id')
    full_df = pd.merge(left=nodes_df, right=node_types_df, how='left', left_on='node_type_id', right_index=True)

    if group_key is not None and len(exclude) > 0:
        # Make sure sure we group-key exists as column
        if group_key not in full_df:
            raise Exception('Could not find column {}'.format(group_key))

        group_keys = set(nodes_df[group_key].unique()) - set(exclude)
        groupings = nodes_df.groupby(group_key)
        # remove any rows with matching column value
        for cond in exclude:
            full_df = full_df[full_df[group_key] != cond]

    return full_df