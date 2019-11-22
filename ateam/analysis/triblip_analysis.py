import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import ateam.analysis.ephys_features.triblip as tri
import pandas as pd


def df_marginal_apply(df, xcol, ycol, data_col, fcn):
    sig_frac_x = df.groupby([xcol])[data_col].apply(fcn)
    sig_frac_x = pd.concat([sig_frac_x], keys=['All'], names=[
                           ycol]).reorder_levels([1, 0])
    sig_frac_y = df.groupby([ycol]).trend_vpost.apply(fcn)
    sig_frac_y = pd.concat([sig_frac_y], keys=['All'], names=[xcol])
    sig_frac = df.groupby([xcol, ycol])[data_col].apply(fcn)
    sig_frac = pd.concat([sig_frac, sig_frac_x, sig_frac_y])
    return sig_frac.unstack(0)


def plot_double_2d(df, xcol, ycol, data_col, fcn_data, fcn_annot, label, **kwargs):
    ax = plt.gca()
    data = df_marginal_apply(df, xcol, ycol, data_col, fcn_data)
    annot = df_marginal_apply(df, xcol, ycol, data_col, fcn_annot)

    sns.heatmap(data, annot=annot, fmt='s', linewidths=0.1, center=0, cmap="RdBu_r",
                cbar_kws={'label': label}, annot_kws={'usetex': False}, **kwargs)
    ax.xaxis.tick_top()
    ax.tick_params('x', labelrotation=45)
    ax.xaxis.label_position = 'top'
    ax.xaxis.labelpad = 20


def plot_sweeps_prop(tri_df, cell, propname, show_burst=True):
    celldat = tri_df.loc[cell]
    celldat = celldat[celldat.is_complete]
    plt.scatter(celldat.frequency, celldat[propname], c='b')
    plt.scatter(celldat[celldat.is_bursty].frequency,
                celldat[celldat.is_bursty][propname], c='r')
    if show_burst:
        plt.legend(['Regular', 'Bursting'])
    for sweep in celldat.index:
        plt.annotate(
            sweep, (celldat.frequency[sweep], celldat[propname][sweep]))
    # plt.ylim(ymin=0)
    plt.ylabel('$V_{post}$ (mV)')
    plt.xlabel('stim frequency (Hz)')


def plot_sweeps_trend(tri_df, cell, show_burst=True):
    return plot_sweeps_prop(tri_df, cell, "postspike_v", show_burst)


def plot_sweeps_aligned_lims(nwb_path, tri_df, cell, dt=0.01):
    for sweepnum in tri_df[tri_df.is_complete == True].loc[cell].dropna().index:
        sweepex = tri.sweepex_from_lims_nwb(nwb_path, sweepnum)
        sweepex.process_spikes()
        t_last = sweepex.spike_feature("peak_t")[-1]
        plt.plot(1000*(sweepex.t - t_last), sweepex.v)
    plt.gca().axvline(1000*dt, color='k')
    plt.xlim(-200, 50)
    plt.xlabel('t (ms)')
    plt.ylabel('V (mV)')


def plot_sweeps_aligned(ctc, tri_df, cell, dt=0.01):
    dataset = ctc.get_ephys_data(cell)
    for sweepnum in tri_df[tri_df.is_complete == True].loc[cell].dropna().index:
        sweep = dataset.get_sweep(sweepnum)
        sweepex = tri.sweepex_from_nwb_sweep(sweep)
        sweepex.process_spikes()
        t_last = sweepex.spike_feature("peak_t")[-1]
        plt.plot(1000*(sweepex.t - t_last), sweepex.v)
    plt.gca().axvline(1000*dt, color='k')
    plt.xlim(-200, 50)
    plt.xlabel('t (ms)')
    plt.ylabel('V (mV)')


def plot_sweep(ctc, cell, sweepnum):
    """Plot a single sweep from web data, using SDK
    """
    dataset = ctc.get_ephys_data(cell)
    sweep = dataset.get_sweep(sweepnum)
    sweepex = tri.sweepex_from_nwb_sweep(sweep)
    sweepex.process_spikes()
    t_last = sweepex.spike_feature("peak_t")[-1]
    t = 1000*(sweepex.t - t_last)

    fig, axes = plt.subplots(2, 1, sharex=True)
    axes[0].plot(t, sweepex.v, color='black')
    axes[1].plot(t, sweepex.i, color='gray')

    axes[0].set_ylabel("mV")
    axes[0].set_title("Voltage")
    axes[1].set_title("Stimulus")
    axes[1].set_ylabel("pA")
    axes[1].set_xlabel("ms")
    plt.xlim(-200, 50)
    # plt.subplots_adjust(hspace=1)


def plot_pop_dist(series_trend, cutoff_z=3):
    bins = np.arange(-5, 15, 0.5)
    series_trend.hist(bins=bins, density=True)
    plt.ylabel('$p(\Delta V_{trend})$', rotation=90)
    plt.xlabel('$\Delta V_{trend}$ (mV)')

    xx = np.linspace(-5, 5, 100)
    popstats = tri.pop_stats_robust(series_trend, cutoff_z=cutoff_z)
    hist_normal = 1/np.sqrt(2*np.pi*popstats.std**2) * \
        np.exp(-(xx-popstats.mean)**2/(2*popstats.std**2))
    plt.plot(xx, hist_normal)
    plt.axvline(popstats.cutoff, color='r')
    plt.legend(['Normal distribution', '{:d}$\sigma$ cutoff'.format(cutoff_z)])
