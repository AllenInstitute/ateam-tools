import pandas as pd
import matplotlib.pyplot as plt

def plot_df_stats_comparison_multiple(df, xvar, yvar, compare, separate):
    df_agg = df.groupby([xvar, compare, separate])[yvar].agg(['mean','std'])
    subsets = df_agg.index.unique(level=separate)
    nplot = len(subsets)
    fig, ax = plt.subplots(nplot)
    for i, name in enumerate(subsets):
        df_sub = df_agg.xs(name, level=separate)
        for name_compare in df_agg.index.unique(level=compare).sort_values():
            df_sub.xs(name_compare, level=compare).plot(y='mean', yerr='std', ax=ax[i], label=name_compare)
        ax[i].set(ylabel=yvar, title=name)

def plot_df_stats_comparison(df, xvar, yvar, compare, ax=None, **plot_args):
    df_agg = df.groupby([xvar, compare])[yvar].agg(['mean','std'])
    ax = ax or plt.axes()
    for name_compare in df_agg.index.unique(level=compare).sort_values():
        df_agg.xs(name_compare, level=compare).plot(y='mean', yerr='std', ax=ax, label=name_compare, **plot_args)
    ax.set(ylabel=yvar)
    return df_agg
