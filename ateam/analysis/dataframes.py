"""Module for analysing and plotting data in pandas dataframes.
"""
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from itertools import combinations

def boxplot(data, var, group, show_swarm=True):
    data = data[~data[var].isna()].sort_values(group)
    nobs = data.groupby(group)[var].count().apply(lambda n: "n={}".format(n))

    ax = sns.boxplot(x=group, y=var, data=data)
    if show_swarm:
        ax2 = sns.swarmplot(x=group, y=var, data=data, color=".01", ax=ax)

    ticks = [tick.get_text() + "\n" + nobs[i] for i, tick in enumerate(ax.get_xticklabels())]
    ax.set_xticklabels(ticks)
    return ax

def boxplot_with_mw_bars(data, var, group, pairs=None, cutoff=0.05, show_swarm=True):
    boxplot(data, var, group, show_swarm=show_swarm)
    pairs, pairs_idx, pvals = pairwise_mw(data, var, group, pairs)
    plot_sig_bars(pvals, pairs_idx, cutoff)

def plot_sig_bars(pvals, pairs_idx, cutoff=0.05):
    ylim = plt.ylim()
    pairs_sig = np.flatnonzero(np.array(pvals)<cutoff)

    yrange = 0.15*(ylim[1]-ylim[0])
    yvals = np.linspace(ylim[1], ylim[1]+yrange, len(pairs_sig))
    for i, pair in enumerate(pairs_sig):
        plot_sig_bar(pvals[pair], yvals[i], pairs_idx[pair])

def plot_sig_bar(pval, y, x_pair):
    plt.plot(x_pair, [y, y], 'k', linewidth=3)
    plt.text(np.mean(x_pair), y, "p={p:.2}".format(p=pval), horizontalalignment='center', verticalalignment='bottom')

def pairwise_mw(data, var, group, pairs=None):
    data = data[~data[var].isna()].sort_values(group)
    group_vals = data[group].unique().tolist()
    group_vals.sort()
    if pairs is None:
        pairs = combinations(group_vals, 2)
    pvals = []
    pairs_idx = []
    groups = data.groupby(group)[var]
    for pair in pairs:
        u, p = mannwhitneyu(groups.get_group(pair[0]), groups.get_group(pair[1]), alternative='two-sided')
        pvals.append(p)
        pairs_idx.append([group_vals.index(pair[0]), group_vals.index(pair[1])])
    return pairs, pairs_idx, pvals