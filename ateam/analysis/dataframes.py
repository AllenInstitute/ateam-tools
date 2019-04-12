"""Module for analysing and plotting data in pandas dataframes.
"""
import warnings
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, linregress
from itertools import combinations
from six import string_types

def flatten_columns(df):
    """Flatten a dataframe with hierarchically indexed columns
    by concatenating labels across levels with underscores 
    """
    df = df.copy()
    df.columns = [col if isinstance(col, string_types) else '_'.join(col).rstrip('_') 
        for col in df.columns.values]
    return df

# TODO: write generic wrapper function?
def group_fit_df(df, xvar, yvar, compare):
#     to use multiple return values in apply(), create Series from dict
    fit_function = lambda df: pd.Series(linfit(*trend_from_df(df, xvar, yvar)))
    df_fit = df.groupby([compare]).apply(fit_function)
    return df_fit

def combine_functions_hierarchical(xvar, yvar, functions, dropna=True):
    # 
    def combined_fcn(df):
        x, y = trend_from_df(df, xvar, yvar, dropna=True)
        outputs = [pd.Series(function(x,y)) for function in functions]
        keys = [function.__name__ for function in functions]
        df_hierarch = pd.concat(outputs, axis=0, keys=keys, names=["analysis", "feature"])
        return df_hierarch
    return combined_fcn

def summary(series):
    unique = series.unique()
    if len(unique)==1:
        out = unique[0]
    else:
        out = "{} values".format(len(unique))
    return out

def trend_from_df(df, xvar, yvar, dropna=True):
    if dropna:
        df = df[~df[yvar].isna()]
    x = df[xvar]
    y = df[yvar]
    return x, y

def threshold(x, y, n_repeats=5):
    error_out = {'min': np.nan, 'mean': np.nan}
    if len(y)<n_repeats:
        return error_out
    x_mins = np.partition(x, n_repeats-1)[:n_repeats]
    return {'min': np.min(x_mins), 'mean': np.mean(x_mins)}
    
def linfit(x, y):
    error_out = {'slope':np.nan, 'yint':np.nan, 'xint':np.nan, 'error':np.nan}
    try:
        if len(y)==0:
            return error_out
        slope, yint = linregress(x,y)[:2]
        xint = np.nan if slope==0 else -yint/slope
        yfit = slope*x + yint
        rmse = np.sqrt(np.mean( (y - yfit) ** 2))
        results = {'slope':slope, 'yint':yint, 'xint':xint, 'error':rmse}
    except Exception as e:
        warnings.warn(e)
        results = error_out
    return results

# Plotting functions
# 

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