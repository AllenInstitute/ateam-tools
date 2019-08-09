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

# Restructuring commands
# 

def flatten_columns(df):
    """Flatten a dataframe with hierarchically indexed columns
    by concatenating labels across levels with underscores 
    """
    df = df.copy()
    df.columns = [col if isinstance(col, string_types) else '_'.join(col).rstrip('_') 
        for col in df.columns.values]
    return df

def summary(series):
    """Aggregation function, summarizes a series by the number of values, or the single shared value if constant
    """
    unique = series.unique()
    if len(unique)==1:
        out = unique[0]
    else:
        out = "{} values".format(len(unique))
    return out


# Plotting functions
# 
def scatterplot_fix(x=None, y=None, data=None, ypad=[], xpad=[], **kwargs):
    """Seaborn scatterplot with limits autoscaling fixed manually
    Default padding 5% of range
    """
    sns.scatterplot(x=x, y=y, data=data, **kwargs)
    plt.xlim(*limits_pad(data[x], *xpad))
    plt.ylim(*limits_pad(data[y], *ypad))

def plot_reg_df(x, y, data=None, **kwargs):
    """Plot regression against a continuous x variable.
    Args use seaborn plotting conventions
    """
    import statsmodels.api as sm
    kwargs.pop('color',None)
    scatterplot_fix(x=x, y=y, data=data, ypad=(0.2, 0), **kwargs)

    xd = data[x]; yd = data[y]
    results = sm.OLS(yd, sm.add_constant(xd), hasconst=True).fit()
    logp = np.log(results.pvalues[1])
    summary = "$R^2={:.2g}$, $log(p)={:.2g}$".format(results.rsquared, logp)
    xpred =  np.linspace(xd.min(), xd.max(), 50)
    ypred = results.predict(sm.add_constant(xpred))
    
    ax = plt.gca()
    ax.plot(xpred, ypred)

    ax.text(0.5, 0.99, summary, transform=ax.transAxes,
            verticalalignment='top', horizontalalignment='center')

def plot_category_df(x, y, data=None, ticks=False, **kwargs):
    """Plot regression against a categorical x variable.
    Args use seaborn plotting conventions
    """
    from statsmodels.formula.api import ols
    results = ols('{} ~ {}'.format(y, x), data=data).fit()
    logp = np.log(results.f_pvalue)
    summary = "$F={:.3g}$, $log(p)={:.3g}$".format(results.fvalue, logp)

    sns.stripplot(x=x, y=y, data=data, palette='muted')
    plt.ylim(*limits_pad(data[y], 0.2, 0))
    ax = plt.gca()
    ax.text(0.5, 0.99, summary, transform=ax.transAxes,
            verticalalignment='top', horizontalalignment='center')
    if not ticks:
        plt.xticks([])

def limits_pad(data, upper=0.05, lower=None):
    xmin = data.min()
    xmax = data.max()
    lower = upper if lower is None else lower
    r = xmax-xmin
    return (xmin - lower*r, xmax + upper*r)

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


# Analysis of two-variable relationships in dataframes
# 
def combine_functions_hierarchical(xvar, yvar, functions, dropna=True):
    """Constructs a single function, to be applied by DataFrame.apply,
    which returns output from several two-variable functions combined into a single hierarchically indexed dataframe.
    """
    def combined_fcn(df):
        x, y = trend_from_df(df, xvar, yvar, dropna=True)
        outputs = [pd.Series(function(x,y)) for function in functions]
        keys = [function.__name__ for function in functions]
        df_hierarch = pd.concat(outputs, axis=0, keys=keys, names=["analysis", "feature"])
        return df_hierarch
    return combined_fcn

def group_fit_df(df, xvar, yvar, compare):
    """Calculate linear fit of two-variable trend independently for each group of data in dataframe 
    
    Arguments:
        xvar, yvar -- columns to fit
        compare -- column to group by
    """
#     to use multiple return values in apply(), create Series from dict
    fit_function = lambda df: pd.Series(linfit(*trend_from_df(df, xvar, yvar)))
    df_fit = df.groupby([compare]).apply(fit_function)
    return df_fit

def trend_from_df(df, xvar, yvar, dropna=True):
    """Extract a pair of columns from dataframe, filtering missing values by default
    """
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