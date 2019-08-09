######################################################
# Authors: Fahimeh Baftizadeh, Taylor Connington
# Date created: 9/1/2017
######################################################

import numpy as np
import pandas as pd
import matplotlib as mlb
import matplotlib.pyplot as plt
import ateam.analysis.brainstim.reobase_utils as ru
from matplotlib.mlab import griddata
import h5py

### EXAMPLES
# SIN
# groupby_cols = ["distance", "fq"]
# agg_cols = ['vm_phase', 'vm_amp', 'vext_phase','vext_amp', 'vi_phase', 'vi_amp']
# merge_cols = ['distance', 'fq']
# t = tpl.get_grouped_dic(cell_list, input_type, stim_type, model_type, inputs, trial, groupby_cols, agg_cols)
# merged = tpl.get_merged_table(t, merge_cols=merge_cols)
# X, Y, Z, Z1, x, y, z, z1 = tpl.get_mesh_X_Y_Z_Z1(merged, xcol="fq", zcol="vm_amp", z1col="vext_amp", ycol="distance")
# ax = tpl.plot_3d_colorbar_mean_mean_scatter(X, Y, Z, Z1, x, y, z, z1, 0, 3)
# ax = tpl.plot_3d_colorbar(X, Y, Z, Z1, 0, 3)

# SIN_DC
# groupby_cols = ["distance", "ic_amp", "fq"]
# agg_cols = ['vm_phase', 'vm_amp', 'vext_phase','vext_amp', 'vi_phase', 'vi_amp', 'vm_stim']
# merge_cols = ['distance', "ic_amp", "fq"]
# t = tpl.get_grouped_dic(cell_list, input_type, stim_type, model_type, inputs, trial, groupby_cols, agg_cols)
# merged = tpl.get_merged_table(t, merge_cols=merge_cols)
# X, Y, Z, Z1, x, y, z, z1 = tpl.get_mesh_X_Y_Z_Z1(merged[merged['fq']==8], xcol="vm_stim", zcol="vm_phase", z1col="vext_phase", ycol="distance")
# ax = tpl.plot_3d_colorbar(X, Y, Z, Z1, 0, 3)

def get_grouped_dic(cell_list, input_type, stim_type, model_type, amp_range, trial, groupby_cols, agg_cols):
    data = {}
    for cell_id in cell_list:
        table = ru.read_cell_tables(input_type=input_type, stim_type=stim_type, model_type=model_type,
                                    cell_gid=cell_id, amp_range=amp_range, trial=trial)
        print "finished reading the table for cell_id:", cell_id
        # index_close_els = ru.get_index_close_els(cell_id, input_type, stim_type, model_type)
        # table = table.drop(index_close_els)
        print "Turned Vm_phase to vm_phase +360"
        table.loc[table["vm_phase"] < 0, "vm_phase"] = table["vm_phase"] + 360
        print "Turned Vext_phase to vext_phase -360"
        table.loc[(table["vext_phase"] < 360) & (table['vext_phase'] > 359), "vext_phase"] = table["vext_phase"] - 360

        data[cell_id] = table.groupby(groupby_cols)[agg_cols].mean().reset_index()
    return data


def get_merged_table(dic, merge_cols):
    dfs = [v for k, v in dic.iteritems()]
    merged = reduce(lambda left, right: pd.merge(left, right, on=merge_cols, how='outer'), dfs)
    return merged


def get_agg_merged_colnames(merged_table, var_name):
    var_agg_merged_colname = np.unique([col for col in merged_table.columns if var_name in col]).tolist()
    return var_agg_merged_colname


def get_stdcol_3d_colorbar(merged_table, var_name):
    col = var_name
    agg_merged_colnames = get_agg_merged_colnames(merged_table, var_name)
    std_col = merged_table[agg_merged_colnames].std(axis=1).as_matrix()
    n_sample = len(std_col)
    std_col = std_col/np.sqrt(n_sample)
    return std_col

def get_meancol_3d_colorbar(merged_table, var_name):
    col = var_name
    agg_merged_colnames = get_agg_merged_colnames(merged_table, var_name)
    mean_col = merged_table[agg_merged_colnames].mean(axis=1).as_matrix()
    return mean_col


def get_mesh_X_Y_Z_Z1(merged_table, xcol, zcol, z1col, ycol="distance", error=False):
    
    x = get_meancol_3d_colorbar(merged_table, xcol)
    z = get_meancol_3d_colorbar(merged_table, zcol)
    z1 = get_meancol_3d_colorbar(merged_table, z1col)
    if error:
        print "error"
        x = get_meancol_3d_colorbar(merged_table, xcol)
        z = get_stdcol_3d_colorbar(merged_table, zcol)
        z1 = get_stdcol_3d_colorbar(merged_table, z1col)

    y = merged_table["distance"].as_matrix()
    
    xi = np.linspace(np.min(x), np.max(x))
    yi = np.linspace(np.min(y), np.max(y))

    X, Y = np.meshgrid(xi, yi)
    Z = griddata(x, y, z, xi, yi, interp='linear')
    Z1 = griddata(x, y, z1, xi, yi, interp='linear')
    return X, Y, Z, Z1, x, y, z, z1


def plot_3d_colorbar(X, Y, Z, Z1, cbar_min, cbar_max, figsize=(18,14), ax=None):
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    import matplotlib.pyplot as plt
    from matplotlib.mlab import griddata
    import numpy as np
    import matplotlib
    import pandas as pd
    fig = plt.figure(figsize=figsize)
    if not ax:
        ax = fig.gca(projection='3d')

    ax.tick_params(axis='x', which='major', pad=5)
    ax.tick_params(axis='y', which='major', pad=5)
    ax.tick_params(axis='z', which='major', pad=20)

    norm = matplotlib.colors.Normalize(vmin=cbar_min,vmax = cbar_max)

    surf1 = ax.plot_surface(X, Y, Z, rstride=1, cstride=1,facecolors=plt.cm.jet(norm(Z)),
                       linewidth=1, antialiased=True)
    surf2 = ax.plot_surface(X, Y, Z1, rstride=1, cstride=1, facecolors=plt.cm.jet(norm(Z1)),
                       linewidth=1, antialiased=True)

    m = cm.ScalarMappable(cmap=plt.cm.jet, norm=norm)
    m.set_array([])
    ax.tick_params(labelsize=35)
    plt.gca().invert_yaxis()
    return ax


def plot_3d_colorbar_mean_mean_scatter(X, Y, Z, Z1, x_mean, y_mean, z_mean, z1_mean, cbar_min, cbar_max, 
                     figsize=(18,14), ax=None):
    '''Same as plot_3d_colorbar except it plots a transparent surface with the points overlayed 
    as a scatter plot. The points represent the mean of mean of the phase or amplitude plotted for
    that distance and frequency/vm stim.'''
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    import matplotlib.pyplot as plt
    from matplotlib.mlab import griddata
    import numpy as np
    import matplotlib
    import pandas as pd
    fig = plt.figure(figsize=figsize)
    if not ax:
        ax = fig.gca(projection='3d')

    ax.tick_params(axis='x', which='major', pad=5)
    ax.tick_params(axis='y', which='major', pad=5)
    ax.tick_params(axis='z', which='major', pad=20)

    norm = matplotlib.colors.Normalize(vmin=cbar_min,vmax = cbar_max)

    # Scatter
    ax.scatter(x_mean, y_mean, z_mean, cmap=plt.cm.jet(norm(Z)), s=50)
    ax.scatter(x_mean, y_mean, z1_mean, cmap=plt.cm.jet(norm(Z)), s=50)

    surf1 = ax.plot_surface(X, Y, Z, rstride=1, cstride=1,facecolors=plt.cm.jet(norm(Z)),
                       linewidth=1, antialiased=True, alpha=0.1)
    surf2 = ax.plot_surface(X, Y, Z1, rstride=1, cstride=1, facecolors=plt.cm.jet(norm(Z1)),
                       linewidth=1, antialiased=True, alpha=0.1)

    m = cm.ScalarMappable(cmap=plt.cm.jet, norm=norm)
    m.set_array([])
    ax.tick_params(labelsize=35)
    plt.gca().invert_yaxis()
    return ax


def plot_mediancol1_col2(gids_list, colname1, colname2, amp, stim_type, model_type, trial):
    med_dic = {}
    d_list = []

    for gid in gids_list:
        t = ru.read_cell_tables(gid, [amp], stim_type, model_type, trial)
        groups = t.groupby([colname2])

        med = []
        for d, g in groups:
            d_list = d_list + [d]
            med.append(g[colname1].agg(np.median, axis=0))

        med_dic[gid] = med

    box_data = [[med_dic[gid][i] for gid in gids_list] for i in range(len(groups))]
    filtered_data = []

    for list in box_data:
        filtered_data.append([x for x in list if ~np.isnan(x)])

    fig = plt.figure()
    fig.set_figheight(7)
    fig.set_figwidth(15)
    ax = fig.add_subplot(111)

    if colname1 == "delta_vm":
        ax.set_title('Median of $\Delta$Vm distribution by distance (amp = {})'.format(amp), fontsize=15)
        ax.set_ylabel('$\Delta$Vm (mV)', fontsize=20)
    else:
        ax.set_title('Median of spike frequency distribution by distance (amp = {})'.format(amp), fontsize=15)
        ax.set_ylabel('Spike frequency (Hz)', fontsize=20)

    ax.set_xlabel('Distance ($\mu$)', fontsize=20)
    boxprops = dict(linewidth=1.5)
    medianprops = dict(linestyle='-', linewidth=1.5, color='firebrick')
    ax.tick_params(labelsize=15)
    ax.boxplot(filtered_data, showmeans=True, boxprops=boxprops, medianprops=medianprops)
    ax.set_xticklabels(["{0:4.2f}".format(x) for x in d_list])
    Fig_folder="/allen/aibs/mat/Fahimehb/Data_cube/reobase/Run_folder/result_figures/dc/progress_report_figs/"
    Fig_name = "Fig.png"
    plt.savefig(Fig_folder + Fig_name)

    plt.show()


def groupby_plot(table, groupby_col, plotx, ploty):
    for groupcol, group in table.groupby(groupby_col):
        fig, ax = plt.subplots()
        ax.plot(group[plotx], group[ploty], marker='o', linestyle='', ms=10)
        ax.set_xlabel(plotx)
        ax.set_ylabel(ploty)
        ax.set_title('{}=  {}'.format(groupby_col, groupcol))
        ax.set_xlim([0,110])
    plt.show()


def passive_comparison(table_pass_peri, table_pass_aa, amp):
    fig, ((ax1, ax2)) = plt.subplots(1, 2, sharex=False, sharey=False)
    fig.set_figheight(5)
    fig.set_figwidth(15)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, hspace=None)

    ax1.plot(table_pass_peri[table_pass_peri["amp"] == amp]["electrode"],
             table_pass_peri[table_pass_peri["amp"] == amp]["delta_deltav"])
    ax2.plot(table_pass_aa[(table_pass_aa["delta_deltav"] <20) & (table_pass_aa["amp"] == amp)]["electrode"],
             table_pass_aa[(table_pass_aa["delta_deltav"] <20) & (table_pass_aa["amp"] == amp)]["delta_deltav"])

    ax1.set_xlabel("electrode_id", size=15)
    ax2.set_xlabel("electrode_id", size=15)
    ax1.set_ylabel("delta_deltav (mV)", size=15)

    ax1.set_title('{}, amp = {}'.format("perisomatic", amp), size=15)
    ax2.set_title('{}, amp = {}'.format("all_active", amp), size=15)
    plt.show()


def Build_subthreshold_comparison_table(gid, Model_Type, Stim_Type, list_trials, inputs):

        table_temp = {}
        table_pass_peri = pd.DataFrame()
        for tr in list_trials:
            table_temp[tr] = ru.read_cell_tables(gid, inputs ,Stim_Type, Model_Type, trial=tr)
            new_deltavm = "delta_vm" + str(tr)
            table_temp[tr] = table_temp[tr].rename(columns={'delta_vm': new_deltavm})
            table_temp[tr]["distance"] = table_temp[tr]["distance"].round()

        table_pass_peri = pd.merge(table_temp[list_trials[0]][["electrode", "x", "y", "z", "distance", "amp", "delta_vm" + str(list_trials[0])]],
                                   table_temp[list_trials[1]][["electrode", "x", "y", "z", "distance", "amp", "delta_vm" + str(list_trials[1]) ]],
                                   on=["electrode", "amp", "x", "y", "z", "distance"]).sort_values(by=["electrode"])
        table_pass_peri["delta_deltav"] = table_pass_peri["delta_vm" + str(list_trials[0])] - table_pass_peri["delta_vm" + str(list_trials[1])]
        # table_pass_peri["avg_deltav"] = (table_pass_peri["delta_vm0"] + table_pass_peri["delta_vm1"]) / 2.
        return table_pass_peri

##################################################
#                                                #
#                 SIN_DC plots                   #
#                                                #
##################################################

def plot_vm_phase_els_Iclamp(table, ic_amp, fq_range):
    fig = plt.figure(figsize=(20,8))
    for fq in fq_range:
        df = table[(table["ic_amp"]==ic_amp) & (table["fq"]==fq)]
        plt.plot(df["electrode"],df["vm_phase"], marker='o', label="vm_fq={}".format(fq))

    plt.title("Vm and Vext phase as a function of el_id for different frequencies at Iclamp={}pA".format(ic_amp*1000), fontsize=20)
    plt.plot(df["electrode"],df["vext_phase"], label="Vext")
    plt.xlabel("Electrode id", fontsize=20)
    plt.ylabel("Phase (degree)", fontsize=20)
    plt.legend(fontsize=15)
    plt.show()


def plot_vm_phase_els_freq(table, ic_amp_range, freq):

    fig = plt.figure(figsize=(20,8))
    for icamp in ic_amp_range:
        df = table[(table["ic_amp"]==icamp) & (table["fq"]==freq)]
        plt.plot(df["electrode"],df["vm_phase"], marker='o', label="vm_icamp={}".format(icamp))

    plt.title("Vm and Vext phase as a function of el_id for different Iclamp at fq={}".format(freq), fontsize=20)
    plt.plot(df["electrode"],df["vext_phase"], label="Vext")
    plt.xlabel("Electrode id", fontsize=20)
    plt.ylabel("Phase (degree)", fontsize=20)
    plt.legend(fontsize=15)
    plt.show()

def plot_vm_amplitude_els_Iclamp(table, ic_amp, fq_range):
    fig = plt.figure(figsize=(20,8))
    for fq in fq_range:
        df = table[(table["ic_amp"]==ic_amp) & (table["fq"]==fq)]
        plt.plot(df["electrode"],df["vm_amp"], marker='o', label="vm_fq={}".format(fq))

    plt.title("Vm and Vext amplitude as a function of el_id for different frequencies at Iclamp={}pA".format(ic_amp*1000), fontsize=20)
    plt.plot(df["electrode"],df["vext_amp"], label="Vext")
    plt.xlabel("Electrode id", fontsize=20)
    plt.ylabel("Amplitude (mV)", fontsize=20)
    plt.legend(fontsize=15)
    plt.show()

def plot_vm_amplitude_els_freq(table, ic_amp_range, freq):

    fig = plt.figure(figsize=(20,8))
    for icamp in ic_amp_range:
        df = table[(table["ic_amp"]==icamp) & (table["fq"]==freq)]
        plt.plot(df["electrode"],df["vm_amp"], marker='o', label="vm_icamp={}".format(icamp))

    plt.title("Vm and Vext amplitude as a function of el_id for different Iclamp at fq={}".format(freq), fontsize=20)
    plt.plot(df["electrode"],df["vext_amp"], label="Vext")
    plt.xlabel("Electrode id", fontsize=20)
    plt.ylabel("Amplitude (mV)", fontsize=20)
    plt.legend(fontsize=15)
    plt.show()


##################################################
#                                                #
#            SIN AND SIN_DC 2D plots             #
#                                                #
##################################################


def get_filtered_table(table, el_dist, icamp_restraint=None, fq_restraint=None):
    # Make sure variables are valid
    if el_dist not in table['distance'].unique():
        raise ValueError('Simulation not performed for electrode at this distance. Should be in microns e.g. 20, 30')
    
    if icamp_restraint is not None:
        if icamp_restraint not in table['ic_amp'].unique():
            raise ValueError('Simulation not performed for electrode at this iclamp value. Should be in nA e.g. 0.03, 0.06 nA')

    if fq_restraint is not None:
        if fq_restraint not in table['fq'].unique():
            raise ValueError('Simulation not performed for electrode at this frequency.')
    
    # Only choose electrodes at a certain distance
    sub_table = table[table['distance'] == el_dist]
    sub_table = sub_table[sub_table['ic_amp'] == icamp_restraint] if icamp_restraint is not None else sub_table
    sub_table = sub_table[sub_table['fq'] == fq_restraint] if fq_restraint is not None else sub_table

    return sub_table


def plot_amp_vs_freq(gid, input_type, stim_type, model_type, amps, trial, el_dist, 
                     icamp=None, ax=None, size=(13,7)):
    if ax is None:
        ax = plt.figure(figsize=size)
        ax = plt.subplot(111)
        
    t = ru.read_cell_tables(gid, amps, input_type, stim_type, model_type, trial)
    sub_table = get_filtered_table(t, el_dist, icamp_restraint=icamp, fq_restraint=None)

    # Get all the stimulation frequencies for x axis
    freqs = sub_table['fq'].unique()
    freqs.sort()
    
    for var_name in ['vext_amp', 'vm_amp', 'vi_amp']:        
        # Initialize an array to store the mean var values (e.g. amp, phase) and err
        data = np.zeros((len(freqs),2))
        
        for i, freq in enumerate(freqs):
            var_data = sub_table[sub_table['fq'] == freq][var_name]
            data[i] = [np.mean(var_data), np.std(var_data)]
        
        # Plot colors to match Costas' paper
        if 'vext' in var_name:
            color = 'r'
        elif 'vm' in var_name:
            color = 'g'
        elif 'vi' in var_name:
            color = 'b'
        
        # Plot
        ax.scatter(freqs, data[:,0], s=25, label=var_name + ', distance = ' + str(el_dist), c=color)
        ax.errorbar(freqs, data[:,0], data[:,1], capsize=6, capthick=3, c=color)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_title('Amplitude vs. frequency for electrode distance {} microns'.format(el_dist), fontsize=14)
    ax.set_xlabel('Frequency', fontsize=14)
    ax.set_ylabel('Amplitude (mV)', fontsize=14)
    ax.legend(loc='center right', fontsize=11)
    
    return ax


def plot_phase_vs_freq(gid, input_type, stim_type, model_type, amps, trial, el_dist, icamp=None, ax=None, size=(13,7)):
    if ax is None:
        ax = plt.figure(figsize=size)
        ax = plt.subplot(111)
    
    t = ru.read_cell_tables(gid, amps, input_type, stim_type, model_type, trial)
    sub_table = get_filtered_table(t, el_dist, icamp_restraint=icamp, fq_restraint=None)
        # Turn phases from 0 to 360
    sub_table.loc[(sub_table["vext_phase"] < 360) & (sub_table['vext_phase'] > 359), "vext_phase"] = sub_table["vext_phase"] - 360
    sub_table.loc[sub_table["vi_phase"] < 0, "vi_phase"] = sub_table["vi_phase"] + 360


    # Get all the stimulation frequencies
    freqs = sub_table['fq'].unique()
    freqs.sort()  
    
    for var_name in ['vext_phase', 'vm_phase', 'vi_phase']: 
        # Initialize an array to store the mean var values (e.g. amp, phase) and err
        data = np.zeros((len(freqs),2))
        
        for i, freq in enumerate(freqs):
            var_data = sub_table[sub_table['fq'] == freq][var_name]
            data[i] = compute_phase_mean_and_std(var_data) 

        # Plot colors to match Costas' paper
        if 'vext' in var_name:
            color = 'r'
        elif 'vm' in var_name:
            color = 'g'
        elif 'vi' in var_name:
            color = 'b'
        
        # Plot
        ax.scatter(freqs, data[:,0], s=25, label=var_name + ', distance = ' + str(el_dist), c=color)
        ax.errorbar(freqs, data[:,0], data[:,1], capsize=6, capthick=3, c=color)

    ax.set_title('Amplitude vs. phase for electrode distance={}'.format(el_dist), fontsize=14)
    ax.set_xlabel('Frequency', fontsize=14)
    ax.set_ylabel('Phase', fontsize=14)
    ax.legend(loc='upper right', fontsize=11)
    ax.set_ylim(-10,370)
    return ax


def compute_phase_mean_and_std(phases):
    '''This function computes the mean and standard deviation for 
    voltage phase data. This is important because the phase is periodic 
    so if there were values around 359 and around 1 it messes up the calculations.'''
    
    def compute_distance_from_0(phase):
        '''e.g. if our phase data is [355,10], we want the avg
        phase to be 15/2=7.5, not 365/2'''
        if phase > 180:
            return 360 - phase
        else:
            return phase
    
    f = np.vectorize(compute_distance_from_0)
    mean = np.mean(f(phases))
    std = np.std(f(phases))
    
    return [mean, std]


def plot_amp_vs_vm_stim(gid, input_type, stim_type, model_type, amps, trial, el_dist, fq, ax=None, size=(13,7)):
    if ax is None:
        ax = plt.figure(figsize=size)
        ax = plt.subplot(111)

    t = ru.read_cell_tables(gid, amps, input_type, stim_type, model_type, trial)
    # Only choose electrodes at a certain distance and frequency
    sub_table = get_filtered_table(t, el_dist, icamp_restraint=None, fq_restraint=fq)
    
    icamps = sub_table['ic_amp'].unique()
    icamps.sort()

    for var_name in ['vext_amp', 'vm_amp', 'vi_amp']:
        data = np.zeros((len(icamps),2))
        vstim = []

        for i, icamp in enumerate(icamps):
            # For each icamp the vstim value is very close to each other, even if the distance changes
            # Assume that the vstim values for each icamp are close enough for plotting
            vstim.append(sub_table[sub_table['ic_amp'] == icamp]['vm_stim'].mean())
            var_data = sub_table[sub_table['ic_amp'] == icamp][var_name]
            data[i] = [np.mean(var_data), np.std(var_data)]

        # Plot colors to match Costas' paper
        if 'vext' in var_name:
            color = 'r'
        elif 'vm' in var_name:
            color = 'g'
        elif 'vi' in var_name:
            color = 'b'
        
        # Plot
        ax.scatter(vstim, data[:,0], s=25, label=var_name + ', distance = ' + str(el_dist), c=color)
        ax.errorbar(vstim, data[:,0], data[:,1], capsize=6, capthick=3, c=color)

    ax.set_title('Amplitude vs. <Vm> for electrode distance={}'.format(el_dist), fontsize=14)
    ax.set_xlabel('<Vm> (mV)', fontsize=14)
    ax.set_ylabel('Amplitude (mV)', fontsize=14)
    ax.legend(loc='center right', fontsize=11)
    return ax


def plot_phase_vs_vm_stim(gid, input_type, stim_type, model_type, amps, trial, el_dist, fq, ax=None, size=(13,7)):
    if ax is None:
        ax = plt.figure(figsize=size)
        ax = plt.subplot(111)

    t = ru.read_cell_tables(gid, amps, input_type, stim_type, model_type, trial)
    # Only choose electrodes at a certain distance and frequency
    sub_table = get_filtered_table(t, el_dist, icamp_restraint=None, fq_restraint=fq)
    # Turn phases from 0 to 360
    sub_table.loc[(sub_table["vext_phase"] < 360) & (sub_table['vext_phase'] > 359), "vext_phase"] = sub_table["vext_phase"] - 360
    sub_table.loc[sub_table["vi_phase"] < 0, "vi_phase"] = sub_table["vi_phase"] + 360

    icamps = sub_table['ic_amp'].unique()
    icamps.sort()

    for var_name in ['vext_phase', 'vm_phase', 'vi_phase']:
        data = np.zeros((len(icamps),2))
        vstim = []

        for i, icamp in enumerate(icamps):
            # For each icamp the vstim value is very close to each other, even if the distance changes
            # Assume that the vstim values for each icamp are close enough for plotting
            vstim.append(sub_table[sub_table['ic_amp'] == icamp]['vm_stim'].mean())
            var_data = sub_table[sub_table['ic_amp'] == icamp][var_name]
            data[i] = compute_phase_mean_and_std(var_data) 
            # data[i] = [np.mean(var_data), np.std(var_data)]

        # Plot colors to match Costas' paper
        if 'vext' in var_name:
            color = 'r'
        elif 'vm' in var_name:
            color = 'g'
        elif 'vi' in var_name:
            color = 'b'
        
        # Plot
        ax.scatter(vstim, data[:,0], s=25, label=var_name + ', distance = ' + str(el_dist), c=color)
        ax.errorbar(vstim, data[:,0], data[:,1], capsize=6, capthick=3, c=color)

    ax.set_title('Phase vs. <Vm> for electrode distance={}'.format(el_dist), fontsize=14)
    ax.set_xlabel('<Vm> (mV)', fontsize=14)
    ax.set_ylabel('Phase', fontsize=14)
    ax.legend(loc='upper right', fontsize=11)
    ax.set_ylim(-10,370)
    return ax


def plot_xcorr(table, input_type, stim_type, model_type, cell_gid, dist, fq, amp, icamp=None, 
               saved_data=None, include_sin=True, ax=None):

    def get_dir_name(el, fq, amp, trial, include_sin, include_iclamp, ic_amp=None):
        """ Find all the outputs for the table"""
        if  not include_sin and not include_iclamp:
            return ru.get_dc_dir_name(el, amp, trial)
        if include_sin and not include_iclamp:
            return ru.get_dir_name(el=el, amp=amp, freq=fq, ic_amp=None, trial=trial)
        if include_iclamp:
            return ru.get_dir_name(el=el, amp=amp, freq=fq, ic_amp=ic_amp, trial=trial)

    
    cell_out_dir = ru.get_output_dir(input_type, stim_type, model_type, cell_gid, saved_data)
    include_iclamp = input_type == ru.InputType.EXTRASTIM_INTRASTIM
    sub_table = plot_helper.get_filtered_table(table, dist, fq_restraint=fq, icamp_restraint=icamp)
    
    sampl_rate = 40000
    maxlags = sampl_rate/fq
            
    ex_delay = 1000
    ex_dur = 9000
    t_start_analysis = ex_delay + 3000.   # 1000 + 3000 = 4000
    t_end_analysis = ex_delay + ex_dur - 2000.   # 8000
    
    ccor_vi = []
    ccor_vm = []
    
    els = sub_table['electrode']
    # print(len(els))
    for el in els:
        dir_name = get_dir_name(el, fq, amp, 0, include_sin, include_iclamp, icamp)
        out_dir =  ru.concat_path(cell_out_dir, dir_name)
        cvh5 = ru.get_cv_files(out_dir, [0])[0]
        dt = cvh5.attrs['dt']

        ve = cvh5['vext'].value
        vm = cvh5['vm'].value
        ve = ve[int(t_start_analysis/dt):int(t_end_analysis/dt)]
        vm = vm[int(t_start_analysis/dt):int(t_end_analysis/dt)]
        vi = np.add(ve, vm)
        
        Nx = len(vi)
        # Compute cross correlation
        c_vi = np.correlate(vi - np.mean(vi), ve - np.mean(ve), mode=2)
        c_vm = np.correlate(vm - np.mean(vm), ve - np.mean(ve), mode=2)

        # Normalize
        c_vi = np.true_divide(c_vi, (Nx * np.std(vi - np.mean(vi)) * np.std(ve - np.mean(ve))))
        c_vm = np.true_divide(c_vm, (Nx * np.std(vm - np.mean(vm)) * np.std(ve - np.mean(ve))))
            
        c_vi = c_vi[Nx - 1 - maxlags:Nx + maxlags]
        c_vm = c_vm[Nx - 1 - maxlags:Nx + maxlags]
        ccor_vi.append(c_vi)
        ccor_vm.append(c_vm)
        
    if ax is None:
        ax = plt.figure(figsize=(11,5))
        ax = plt.subplot(111)
    t = np.linspace(-1./fq, 1./fq, len(ccor_vi[0]))
    ax.plot(t, np.mean(ccor_vi, axis=0), c='b', label='vi_ve')
    ax.plot(t, np.mean(ccor_vm, axis=0), c='g', label='vm_ve')
    ax.errorbar(t, np.mean(ccor_vi, axis=0), np.std(ccor_vi, axis=0), alpha=0.01, c='b')
    ax.errorbar(t, np.mean(ccor_vm, axis=0), np.std(ccor_vm, axis=0), alpha=0.01, c='g')    
    ax.set_xlim(-1./fq, 1./fq)
    ax.set_xticks(np.linspace(-1./fq, 1./fq, 3))

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
#     ax.set_xlabel('Time period (s)')
#     ax.set_ylabel('xcorr')


##################################################
#                                                #
#                 3D plots                       #
#                                                #
##################################################

# groupby_cols = ['distance', 'fq']
# var_cols = ['vm_amp', 'vm_phase', 'vext_amp', 'vext_phase', 'vi_amp', 'vi_phase']

# groupby_cols_sindc = ['distance', 'fq', 'ic_amp']
# var_cols_sindc = ['vm_amp', 'vm_phase', 'vext_amp', 'vext_phase', 'vi_amp', 'vi_phase', 'vm_stim']

### EXAMPLE
# t = get_concat_df(cell_list, 'extrastim', 'sin', 'all_active', [-0.0002], 0)
# t_merged = get_merged_mean(t, groupby_cols, var_cols)
# plot_3d_colorbar_transparent(t_merged, 'fq', 'vm_amp', 'vext_amp', 0, 3)

# t = get_concat_df(cell_list, 'extrastim_intrastim', 'sin_dc', 'all_active', [-0.0002], 0)
# t_merged = get_merged_mean(t, groupby_cols_sindc, var_cols_sindc)
# plot_3d_colorbar_transparent(t_merged[t_merged['fq']==8], 'vm_stim', 'vm_amp', 'vext_amp', 0, 3)


def get_concat_df(cell_list, input_type, stim_type, model_type, amp_range, trial):
    df = pd.DataFrame()
    for cell_id in cell_list:
        table = ru.read_cell_tables(input_type=input_type, stim_type=stim_type, model_type=model_type,
                                    cell_gid=cell_id, amp_range=amp_range, trial=trial)
        print "finished reading the table for cell_id:", cell_id
        # Move vext phase from 360 to 0
        table.loc[(table["vext_phase"] < 360) & (table['vext_phase'] > 359), "vext_phase"] = table["vext_phase"] - 360
        df = pd.concat([df, table])
    return df


def get_merged_mean(df, groupby_cols, var_cols):
    return df.groupby(groupby_cols)[var_cols].mean().reset_index()


def get_merged_sem(df, groupby_cols, var_cols):
    return df.groupby(groupby_cols)[var_cols].sem().reset_index()


def get_mesh(table, xcol, zcol, z1col, ycol='distance'):
    x = table[xcol].values
    z = table[zcol].values
    z1 = table[z1col].values
    y = table[ycol].values
    
    xi = np.linspace(np.min(x), np.max(x))
    yi = np.linspace(np.min(y), np.max(y))
    
    X, Y = np.meshgrid(xi, yi)
    Z = griddata(x, y, z, xi, yi, interp='linear')
    Z1 = griddata(x, y, z1, xi, yi, interp='linear')
    return X, Y, Z, Z1


def plot_3d_colorbar_transparent(merged_table, xcol, zcol, zcol1, cbar_min, cbar_max, ycol='distance', size=(18,10)):
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    import matplotlib.pyplot as plt
    from matplotlib.mlab import griddata
    import matplotlib
    import numpy as np
    import pandas as pd
    
    fig = plt.figure(figsize=size)
    ax = fig.gca(projection='3d')
    
    X, Y, Z, Z1 = get_mesh(merged_table, xcol, zcol, zcol1)
    norm = matplotlib.colors.Normalize(vmin=cbar_min, vmax=cbar_max)
    
    ax.tick_params(axis='x', which='major', pad=5)
    ax.tick_params(axis='y', which='major', pad=5)
    ax.tick_params(axis='z', which='major', pad=20)

    # Scatter
    ax.scatter(merged_table[xcol], merged_table[ycol], merged_table[zcol], cmap=plt.cm.jet(norm(Z)), s=50)
    ax.scatter(merged_table[xcol], merged_table[ycol], merged_table[zcol1], cmap=plt.cm.jet(norm(Z)), s=50)
    
    # Surface
    surf1 = ax.plot_surface(X, Y, Z, rstride=1, cstride=1,facecolors=plt.cm.jet(norm(Z)),
                            linewidth=1, antialiased=True, alpha=0.1)
    surf2 = ax.plot_surface(X, Y, Z1, rstride=1, cstride=1, facecolors=plt.cm.jet(norm(Z1)),
                            linewidth=1, antialiased=True, alpha=0.1)
    
    m = cm.ScalarMappable(cmap=plt.cm.jet, norm=norm)
    m.set_array([])
    #cbar=plt.colorbar(m)
    #cbar.ax.tick_params(labelsize=20)
    ax.tick_params(labelsize=35)
    plt.gca().invert_yaxis()
    return ax


# For sin_dc phase plot
def plot_3d_colorbar_transparent2(merged_table, xcol, zcol, zcol1, cbar_min, cbar_max, ycol='distance', size=(18,10)):
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    import matplotlib.pyplot as plt
    from matplotlib.mlab import griddata
    import matplotlib
    import numpy as np
    import pandas as pd
    
    fig = plt.figure(figsize=size)
    ax = fig.gca(projection='3d')
    
    X, Y, Z, Z1 = get_mesh(merged_table, xcol, zcol, zcol1)
    norm = matplotlib.colors.Normalize(vmin=cbar_min, vmax=cbar_max)
    
    ax.tick_params(axis='x', which='major', pad=5)
    ax.tick_params(axis='y', which='major', pad=5)
    ax.tick_params(axis='z', which='major', pad=20)

    # Scatter
    ax.scatter(merged_table[xcol], merged_table[ycol], merged_table[zcol], cmap=plt.cm.jet(norm(Z)), s=50)
    ax.scatter(merged_table[xcol], merged_table[ycol], merged_table[zcol1], cmap=plt.cm.jet(norm(Z)), s=50)
    
    # Surface
    surf2 = ax.plot_surface(X, Y, Z1, rstride=1, cstride=1, facecolors=plt.cm.jet(norm(Z1)),
                            linewidth=1, antialiased=True, alpha=0.1)
    surf1 = ax.plot_surface(X, Y, Z, rstride=1, cstride=1,facecolors=plt.cm.jet(norm(Z)),
                            linewidth=1, antialiased=True, alpha=0.1)
    
    m = cm.ScalarMappable(cmap=plt.cm.jet, norm=norm)
    m.set_array([])
    #cbar=plt.colorbar(m)
    #cbar.ax.tick_params(labelsize=20)
    ax.tick_params(labelsize=35)
    plt.gca().invert_yaxis()
    return ax


################################################

#                  FFT                         #

#################################################

def plot_fft(sig_fft, sample_freq , ax=None):
    if ax is None:
        ax = plt.figure(figsize=(20,5))
        ax = plt.subplot(111)

    sorted_sig_fft = np.asarray([x for _,x in sorted(zip(sample_freq, sig_fft))])
    sorted_sampling_freq = np.sort(sample_freq)
    power = np.abs(sorted_sig_fft)
    ax.semilogy(sorted_sampling_freq, np.abs(power), label='fft')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('semilog power')
    ax.set_xlim(-110,110)
    ax.legend(loc='best')
    return ax


################################################

#                  NWB                         #

#################################################

def plot_nwb_trace(v_trace, sampling_freq, title=None, ax=None, label=None):
    if ax is None:
        ax = plt.figure(figsize=(20,5))
        ax = plt.subplot(111)

    N = len(v_trace)
    dt = 1. / (sampling_freq * 1000)
    tstop = N * dt
    time = np.arange(0, tstop, dt)
    ax.plot(time, v_trace, label=label)
    ax.set_xlabel('Time(s)', size=15)
    ax.set_ylabel('Voltage(mV)', size =15)
    ax.legend(loc='best', fontsize=15)
    if title:
        ax.set_title(title, size= 20)
    return ax


def get_cut_window(in_amp, ex_delay, ex_dur):
    if in_amp == 0:
        t_start = ex_delay + 500
        t_end = ex_delay + ex_dur - 1000
    else:
        t_start = ex_delay + 2500
        t_end = ex_delay + ex_dur - 2000
    return t_start, t_end

def filter_list(in_amp, ex_delay, ex_dur,var1_list, var2_list):
    t_start , t_end = get_cut_window(in_amp, ex_delay, ex_dur)
    ndx = np.where((var1_list >= t_start) & (var1_list <= t_end))
    l = [var2_list[i] for i in ndx]
    flat_list = [item for sublist in l for item in sublist]
    return flat_list

def phase_correction(spike_phase_A_list):
        temp = [x + (1.5 * np.pi) for x in spike_phase_A_list]
        return [(x / (2 * np.pi) - int(x / (2 * np.pi))) * 2 * np.pi for x in temp]

# def filter_list(row, **kwargs):
#     in_amp = row['in_amp(pA)']
#     ex_delay = row['ex_delay(ms)']
#     ex_dur = row['ex_dur(ms)']
#     t_start , t_end = get_cut_window(in_amp, ex_delay, ex_dur)
#     varname1 = kwargs['varname1']
#     varname2 = kwargs['varname2']
#     ndx = np.where((row[varname1] >= t_start) & (row[varname1] <= t_end))
#     l = [row[varname2][i] for i in ndx]
#     flat_list = [item for sublist in l for item in sublist]
#     return flat_list

# def phase_correction(row):
#     temp = [x + (1.5 * np.pi) for x in row['spike_phase_A']]
#     return [(x/(2*np.pi) - int(x/(2*np.pi))) * 2 * np.pi for x in temp]