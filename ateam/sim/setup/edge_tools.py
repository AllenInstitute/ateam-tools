import random, numpy as np

def rand_connector(source, target, prob, nsyn_min=0, nsyn_max=10):
    """
    Connect source and target pair with given probability
    returns random number of synapses chosen from range.
    For use with default iterator (one-to-one)
    """
    if source.node_id == target.node_id or random.random() > prob:
        return None
    tmp_nsyn = random.randint(nsyn_min, nsyn_max)
    return tmp_nsyn


def select_source_cells(sources, target, nsources_min, nsources_max, nsyn_min=0, nsyn_max=10):
    """
    Connect target to a random subset of sources.
    For use with all-to-one iterator.
    """
    total_sources = len(sources)
    nsources = np.random.randint(nsources_min, nsources_max)
    selected_sources = np.random.choice(total_sources, nsources, replace=False)
    syns = np.zeros(total_sources)
    syns[selected_sources] = np.random.randint(nsyn_min, nsyn_max, size=nsources)
    return syns