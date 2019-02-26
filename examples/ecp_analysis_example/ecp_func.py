import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import icsd
import neo
import quantities as pq
import csv
import pylab as pl
from scipy.interpolate import spline
from allensdk.core.cell_types_cache import CellTypesCache
import math
from scipy import signal

def get_loc_electrodes(filepath,filename):
	#"""
	#Read the location of electrode
	#"""
	with open(filepath+filename) as csvDataFile:
	    csvReader = csv.reader(csvDataFile, delimiter=' ')
	    headers=csvReader.next()
	    #for row in csvReader:
	    #    print(row)
	    
	    column = {}
	    for h in headers:
		column[h] = []
		
	    for row in csvReader:
		for h, v in zip(headers, row):
		    column[h].append(float(v))
		    
	z=column['z_pos']   
	x=column['x_pos']
	y=column['y_pos']
        return (x,y,z)

# calculate the spread of amplitude 
def cal_spread(A0,dist0,th=0.12):    
    dist = np.linspace(dist0.min(),dist0.max(),10000) #1000 represents number of points to make between T.min and T.max
    A = spline(dist0,A0,dist)

    midamp = max(A)*th
    indx = range(0,A.argmax())
    idx1 = find_nearest(A[indx],midamp)
    indx = range(A.argmax(),len(A))
    idx2 = find_nearest(A[indx],midamp)+A.argmax()
    return dist[idx2]-dist[idx1] 


def simpleaxis(ax):
    #Hide the right and top spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx    

def find_Onset_time(waveform0,t_STA0):
    t_STA = np.linspace(t_STA0.min(),t_STA0.max(),10000) #300 represents number of points to make between T.min and T.max
    waveform = np.array(spline(t_STA0,waveform0,t_STA))

    rate=np.array([0.0]*(len(waveform))) 
    for i in range(0,len(waveform)-1):
        rate[i]= (waveform[i+1]-waveform[i])/dt/1000
    return t_STA[rate.argmin()]

def find_neg_peak_time(waveform0,t_STA0):
    t_STA = np.linspace(t_STA0.min(),t_STA0.max(),10000) #300 represents number of points to make between T.min and T.max
    waveform = np.array(spline(t_STA0,waveform0,t_STA)) 
    return t_STA[find_neg_peak_indx(waveform,t_STA)]

def find_neg_peak_indx(waveform,t_STA):
    indx = np.where(np.logical_and(t_STA>=-0.2,t_STA<=1.5))
    t_STA= t_STA[indx]
    waveform = waveform[indx]
    return waveform.argmin()+indx[0][0]

# trough to peak width
def calc_tp_width(waveform0,t_STA0):
    
    t_STA = np.linspace(t_STA0.min(),t_STA0.max(),10000) #10000 represents number of points to make between T.min and T.max
    waveform = spline(t_STA0,waveform0,t_STA)

    idx1 = find_neg_peak_indx(waveform,t_STA) 

    indx = range(idx1,len(waveform))
    idx2 = waveform[indx].argmax()+idx1
    
    return t_STA[idx2]-t_STA[idx1]

#half-way width
def calc_width(waveform0,t_STA0):
    
    t_STA = np.linspace(t_STA0.min(),t_STA0.max(),10000) #10000 represents number of points to make between T.min and T.max
    waveform = spline(t_STA0,waveform0,t_STA)
    
    NegIndx = find_neg_peak_indx(waveform,t_STA) 
    
    indx = range(0,NegIndx)
    firstPeakIndx = waveform[indx].argmax()
    indx = range(NegIndx,len(waveform))
    secondPeakIndx = waveform[indx].argmax()+NegIndx

    if secondPeakIndx==NegIndx:
        secondPeakIndx=len(waveform)
    
    midamp = (waveform[0]-waveform[NegIndx])*0.5+waveform[NegIndx]
   
    indx1 = range(firstPeakIndx,NegIndx)
    idx1 = find_nearest(waveform[indx1],midamp)+firstPeakIndx
    
    indx2 = range(NegIndx,secondPeakIndx)
    idx2 = find_nearest(waveform[indx2],midamp)+NegIndx
    
    return t_STA[idx2]-t_STA[idx1]


def calc_cap(waveform):
    indx = range(0,waveform.argmin())
    cap = max(waveform[indx])-waveform[0]
    return cap


def get_icsd(data,fs,dh,method='spline'): 
    # data: row--channels,column--time series; fs:sampling frequency;  dh: the distance between each channels
    #prepare lfp data for use, by changing the units to SI and append quantities,
    #along with electrode geometry, conductivities and assumed source geometry
    nch=data.shape[0] #the number of channels
    lfp_data = data * 1E-6 * pq.V        # [uV] -> [V]
    z_data = np.linspace(0E-6, dh*nch*1E-6, nch) * pq.m 
    diam = 40E-6 * pq.m                                  # [m]
    h = dh*1E-6 * pq.m                                   # [m]
    sigma = 0.3 * pq.S / pq.m                            # [S/m] or [1/(ohm*m)]

    lfp_data = neo.AnalogSignal(lfp_data.T, sampling_rate=fs/1000*pq.kHz)
    # Input dictionaries for each method
    delta_input = {
        'lfp' : lfp_data,
        'coord_electrode' : z_data,
        'diam' : diam,          # source diameter
        'sigma' : sigma,        # extracellular conductivity
        'sigma_top' : sigma,    # conductivity on top of cortex
        'f_type' : 'gaussian',  # gaussian filter
        'f_order' : (3, 1),     # 3-point filter, sigma = 1.
    }
    step_input = {
        'lfp' : lfp_data,
        'coord_electrode' : z_data,
        'diam' : diam,
        'h' : h,                # source thickness
        'sigma' : sigma,
        'sigma_top' : sigma,
        'tol' : 1E-12,          # Tolerance in numerical integration
        'f_type' : 'gaussian',
        'f_order' : (3, 1),
    }
    spline_input = {
        'lfp' : lfp_data,
        'coord_electrode' : z_data,
        'diam' : diam,
        'sigma' : sigma,
        'sigma_top' : sigma,
        'num_steps' : 201,      # Spatial CSD upsampling to N steps
        'tol' : 1E-12,
        'f_type' : 'gaussian',
        'f_order' : (20, 5),
    }

    std_input = {
        'lfp' : lfp_data,
        'coord_electrode' : z_data,
        'sigma' : sigma,
        'f_type' : 'gaussian',
        'f_order' : (3, 1),
    }

    if method=='delta':
        csd_dict = icsd.estimate_csd(**delta_input)
    elif method=='step':
        csd_dict = icsd.estimate_csd(**step_input)
    elif method=='std':
        csd_dict = icsd.estimate_csd(**std_input)
    else: #default method =='spline'
        csd_dict = icsd.estimate_csd(**spline_input)
    
    csd = csd_dict[0].magnitude.T   #original
    csd_smooth = csd_dict[1].magnitude.T  #filtered

    return (csd,csd_smooth)

def butter_highpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_highpass_filter(data, cutoff, fs, order=5):
    b, a = butter_highpass(cutoff, fs, order=order)
    y = signal.filtfilt(b, a, data)
    return y

def butter_lowpass(cutoff, fs, order=3):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=3):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = signal.filtfilt(b, a, data)
    return y

def plot_rec_data(samplingrate,lgn_spikes,Vi,spikes,Ve,fig_size=(3.5,5.5)):
    
    T=Vi.shape[0] #The total data points in time serires
    dt = 1000.0/samplingrate # from m to ms
    t=pl.frange(0,T*dt-dt,dt)

    # plot with various axes scales
    fig, ax = plt.subplots(figsize=fig_size)

    # Plot lgn input
    ax=plt.subplot(411)
    plt.eventplot(lgn_spikes,color='k') 
    plt.ylabel('LGN input')
    plt.xlim([0,T*dt])    
    simpleaxis(ax)
    
    # Plot membrane potential
    ax=plt.subplot(412)
    plt.plot(t,Vi,color='k')
    plt.ylabel('Vi (mV)')
    locs, labels = plt.yticks()
    ax.set(yticks=np.arange(locs[0],locs[-1]+1,(locs[-1]-locs[0])/5))
    simpleaxis(ax)
    
    # plot spikes
    ax=plt.subplot(413)
    firingrate = len(spikes)/(round(T*dt))*1000
    ax.eventplot(spikes,color='k')
    plt.ylabel('spikes')
    plt.yticks([])
    plt.xlim([0,T*dt])
    plt.title('Firing rate='+str(np.around(firingrate,decimals=2))+'Hz')
    simpleaxis(ax)

    
    #plot ecp
    ax=plt.subplot(414)
    somaindx = np.unravel_index(Ve.argmin(),Ve.shape)  #The recording sites that have the largest negative amplitude
    plt.plot(t,Ve[somaindx[0],:],color='k')
    plt.ylabel('Ve (${\mu}$V)')
    simpleaxis(ax)
    plt.xlim([0,T*dt])
    plt.xlabel('Time (ms)')
    locs, labels = plt.yticks()
    ax.set(yticks=np.arange(locs[0],locs[-1]+1,(locs[-1]-locs[0])/4))

    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.6,wspace=0.3)

    return fig

# function to calculate spike triggered average and covariance for one channels
def cal_STA_STC_onech(data,spikes,win,samplingrate): # data: one_dimensional data
    lwin = win[0]
    hwin = win[1]
    dt = 1000.0/samplingrate # from m to ms
    t_STA = pl.frange(lwin,hwin-dt,dt)
    
    for i in range(0,len(spikes)):
        indx = range(int(spikes[i,]/dt)+int(lwin/dt),int(spikes[i,]/dt)+int(hwin/dt))   #-1ms, +3ms
        if i==0:
            arr=np.array(data[indx])
        else:
            if len(data)>int(spikes[i,]/dt)+int(hwin/dt):
                arr = np.vstack((arr, np.array(data[indx])))
    return t_STA,np.mean(arr,axis=0),np.std(arr,axis=0)

# function to calculate spike triggered average and covariance for multiple channels
def cal_STA_STC(data,spikes,win,samplingrate): # ecp_data: channels*times; output: channels*ecp_win_times
    
    num_channels = data.shape[0]  # number of channels
    
    dt = 1000.0/samplingrate      # from m to ms
    times=int((win[1]-win[0])/dt)
    
    #Calculate STA and STC for multiple channels
    STA = np.array([[0.0]*times]*num_channels)     #spike triggered average
    STC = np.array([[0.0]*times]*num_channels)     #spike triggered standard deviation 

    for i in range(num_channels):
        t_STA,STA[i,:],STC[i,:]=cal_STA_STC_onech(data[i,:],spikes,win,samplingrate)
    
    return (t_STA,STA,STC)

def cal_EAP_features(t_STA,Ve_STA,Th):
    num_channels = Ve_STA.shape[0]
    A = np.zeros(num_channels)
    W =  np.zeros(num_channels)
    TPW =  np.zeros(num_channels)
    t_NegPeak =  np.zeros(num_channels)
    for i in range(num_channels): 
        ecp_STA = Ve_STA[i,:]
        ecp_STA = ecp_STA-ecp_STA[0]   
        if abs(min(ecp_STA))>=Th:
            A[i] = abs(min(ecp_STA))    
            W[i] = calc_width(ecp_STA,t_STA)
            TPW[i] = calc_tp_width(ecp_STA,t_STA)
            t_NegPeak[i]= find_neg_peak_time(ecp_STA,t_STA)  #negative peak time
    return A,W,TPW,t_NegPeak

def data_reshape(data,eleRowN,eleColN):
    output = np.array([[0.0]*(eleColN)]*(eleRowN))
    for jj in range(eleColN): 
        for ii in range(eleRowN):
            output[ii,jj] = data[ii+eleRowN*jj]
    return output

def plot2D_features(morph,xx,yy,zA_mask,zW_mask,zTPW_mask,ztNP_mask,fig_size=(10.5,1.6)):
    fig, ax = plt.subplots(figsize=fig_size)
    ax=plt.subplot(141) 
    plot_cell_morph_xy(ax,morph,1)
    heatmap = plt.pcolor(xx,yy,zA_mask,zorder=20,alpha=0.7,edgecolor=(1,1,1,0),linewidth=0,rasterized=True)
    plt.colorbar()
    plt.title('Amplitude (${\mu}$V)')
    ax.set_xlim([-80,80])
    ax.set_ylim([-530,-370])
    clean_axis(ax)
    
    ax=plt.subplot(142)
    plot_cell_morph_xy(ax,morph,1)
    heatmap = plt.pcolor(xx,yy,zW_mask,zorder=20,alpha=0.7,edgecolor=(1,1,1,0),linewidth=0,rasterized=True)#,vmin=0.2, vmax=0.7)
    plt.colorbar()
    ax.set_xlim([-80,80])
    ax.set_ylim([-530,-370])
    plt.title('Width (ms)')
    clean_axis(ax)
    
    
    ax=plt.subplot(143)
    plot_cell_morph_xy(ax,morph,1)
    heatmap = plt.pcolor(xx,yy,zTPW_mask,zorder=20,alpha=0.7,edgecolor=(1,1,1,0),linewidth=0,rasterized=True)#,vmin=0.2, vmax=0.7)
    plt.colorbar()
    ax.set_xlim([-80,80])
    ax.set_ylim([-530,-370])
    plt.title('TP width (ms)')
    clean_axis(ax)
    
    ax=plt.subplot(144)
    plot_cell_morph_xy(ax,morph,1)
    heatmap = plt.pcolor(xx,yy,ztNP_mask,zorder=20,alpha=0.7,edgecolor=(1,1,1,0),linewidth=0,rasterized=True)#,vmin=-0.2, vmax=0.6)
    plt.colorbar()
    ax.set_xlim([-80,80])
    ax.set_ylim([-530,-370])
    plt.title('Peak latency (ms)')
    clean_axis(ax)
    
    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.5,
                wspace=0.5)
    
    return fig




def plot_cell_morph_xy(axes,morph,black=0):
    
    if black==1:
        soma_col='grey'#'crimson'#[194.0/255.0,67.0/255.0,55.0/255.0]
        axon_col='dimgrey'
        dend_col='dimgrey'
        apical_dend_col='k'
        lw=0.5
    else:
        soma_col='crimson'#[194.0/255.0,67.0/255.0,55.0/255.0]
        axon_col=[93.0/255.0,127.0/255.0,177.0/255.0]
        dend_col=[153.0/255.0,40.0/255.0,39.0/255.0]
        apical_dend_col=[227.0/255.0,126.0/255.0,39.0/255.0]
        lw=1#0.5
    
    for n in morph.compartment_list:
        for c in morph.children_of(n):
            if n['type']==2:
                axes.plot([n['x'], c['x']], [n['y'], c['y']], color=axon_col,linewidth=lw)

            if n['type']==3:
                axes.plot([n['x'], c['x']], [n['y'], c['y']], color=dend_col,linewidth=lw)
     
            if n['type']==4:
                axes.plot([n['x'], c['x']], [n['y'], c['y']], color=apical_dend_col,linewidth=lw)
           
            if n['type']==1: #soma
                axes.scatter(n['x'],n['y'],s=math.pi*(n['radius']**2),color=soma_col,zorder=3)
   # axes.set_ylabel('y')
   # axes.set_xlabel('x')
    simpleaxis(axes)
    
    
# function to plot morphology
def cell_morphology_rot(cell_id, x_soma, y_soma, z_soma, theta):

    theta_z = theta[2]
    theta_y = theta[1]
    theta_x = theta[0]
    
    # download and open an SWC file
    ctc = CellTypesCache(manifest_file='cell_types/manifest.json')
    morph = ctc.get_reconstruction(cell_id) 
    
    #First applying a rotation angle around z axis
    tr_rot_z = [math.cos(theta_z),-math.sin(theta_z),0,
                math.sin(theta_z),math.cos(theta_z),0,
                0,0,1,
                0,0,0
              ]   
    
    #Second applying a rotation angle around y axis
    tr_rot_y = [math.cos(theta_y),0,math.sin(theta_y),
                0,1,0,
                -math.sin(theta_y),0,math.cos(theta_y),
                0,0,0
              ]
    
    #Third applying a rotation angle around x axis
    tr_rot_x = [1,0,0,
                0,math.cos(theta_x),-math.sin(theta_x),            
                0,math.sin(theta_x),math.cos(theta_x),
                0,0,0
              ]

    
    morph.apply_affine(tr_rot_x)
    morph.apply_affine(tr_rot_y)
    morph.apply_affine(tr_rot_z)
       
    # translate the soma location
    tr_soma = [1, 0, 0,
               0, 1, 0,
               0, 0, 1,
               -morph.soma["x"]+x_soma, -morph.soma["y"]+y_soma, -morph.soma["z"]+z_soma
             ]
    morph.apply_affine(tr_soma)
    

    # Make a line drawing of x-y and y-z views    
    return morph

    

def plot_EAP(t_STA,STA,fig_size=(5.5,3.5)):  
    
    somaindx = np.unravel_index(STA.argmin(),STA.shape)  #The recording sites that have the largest negative amplitude
    maxCh=somaindx[0]
    fig, ax = plt.subplots(figsize=fig_size)
    ax =plt.subplot(121)
    for i in range(maxCh-5,maxCh+6): 
        if i==maxCh:
            plt.plot(t_STA,STA[i,:]*0.01+i,color = 'grey',linewidth = 3 )  
        else:
            plt.plot(t_STA,STA[i,:]*0.01+i,color='grey',linewidth = 1.5 )  
    plt.ylim([maxCh-7,maxCh+7])
    plt.xlabel('Time (ms)')   
    plt.title('EAP') 
    plt.yticks([])
    plt.xlim([t_STA[0],t_STA[-1]+0.1])
    ax.set(xticks=np.arange(t_STA[0],t_STA[-1]+0.1,1))
    simpleaxis(ax)


    ax=plt.subplot(122)
    for i in range(maxCh-5,maxCh+6): 
        if i==maxCh:
            plt.plot(t_STA,STA[i,:]/abs(min(STA[i,:]))+i,color = 'grey',linewidth = 3 )  
        else:
            plt.plot(t_STA,STA[i,:]/abs(min(STA[i,:]))+i,color='grey',linewidth = 1.5 )  
    plt.ylim([maxCh-7,maxCh+7])
    plt.xlabel('Time (ms)')   
    plt.title('Normalized EAP') 
    plt.yticks([])
    plt.xlim([t_STA[0],t_STA[-1]+0.1])
    ax.set(xticks=np.arange(t_STA[0],t_STA[-1]+0.1,1))
    simpleaxis(ax)

    return fig



def plot1D_features(dist,A,W,TPW,t_NegPeak,maxd=50,fig_size=(10.5,1.5)):  
    
    maxCh=np.unravel_index(A.argmax(),A.shape)[0] #The recording sites that have the largest negative amplitude
    
    fig, ax = plt.subplots(figsize=fig_size)
  
    ax=plt.subplot(141)
    plt.plot(dist,A,'o-',color='k')
    plt.ylabel('Amplitude (${\mu}$V)')
    plt.xlabel('Distance (${\mu}$m)')
    simpleaxis(ax)
    plt.xlim([-maxd,maxd+1])
    ax.set(xticks=np.arange(-maxd,maxd+1,maxd))
    locs, labels = plt.yticks()
    ax.set(yticks=np.arange(locs[0],locs[-1]+5,(locs[-1]-locs[0])/2))
    
    ax=plt.subplot(142)   
    plt.plot(dist,W,'o-',color = 'k')
    plt.ylabel('Width (ms)')
    plt.xlabel('Distance (${\mu}$m)')
    simpleaxis(ax)
    locs, labels = plt.yticks()
    ax.set(yticks=np.arange(locs[0],locs[-1]+0.1,(locs[-1]-locs[0])/2))
    plt.xlim([-maxd,maxd+1])
    ax.set(xticks=np.arange(-maxd,maxd+1,maxd))

    ax=plt.subplot(143)   
    plt.plot(dist,TPW,'o-',color = 'k')
    plt.ylabel('TP Width (ms)')
    plt.xlabel('Distance (${\mu}$m)')
    simpleaxis(ax)
    locs, labels = plt.yticks()
    ax.set(yticks=np.arange(locs[0],locs[-1]+0.1,(locs[-1]-locs[0])/2))
    plt.xlim([-maxd,maxd+1])
    ax.set(xticks=np.arange(-maxd,maxd+1,maxd))
    
    
    ax=plt.subplot(144) 
    plt.plot(dist,t_NegPeak,'o-',color = 'k')
    plt.ylabel('Peak latency (ms)')    
    plt.xlabel('Distance (${\mu}$m)')
    simpleaxis(ax)
    locs, labels = plt.yticks()
    plt.ylim([locs[0],locs[-1]+0.1])
    ax.set(yticks=np.arange(locs[0],locs[-1]+0.1,(locs[-1]-locs[0])/2))
    plt.xlim([-maxd,maxd+1])
    ax.set(xticks=np.arange(-maxd,maxd+1,maxd))
    
    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.6,
                    wspace=0.6)
    return fig

def clean_axis(ax):
    #Hide the right and top spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()       
    ax.set_yticks([])
    ax.set_xticks([])
    ax.spines['bottom'].set_visible(False)

def plot_lfp_csd(lfp,csd,xrange,yrange,fs,fig_size=(5.5,3.5)):
    indx1 = xrange[0]*fs
    indx2 = xrange[1]*fs

    fig, ax = plt.subplots(figsize=fig_size)  
    ax=plt.subplot(121)  
    for ch in range(yrange[0],yrange[1]):
        chunkch = lfp[ch,:]
        plt.plot(chunkch+ch*100,'k',lw=.5)
    plt.xlabel('Time (ms)')  
    plt.xlim([indx1,indx2+1])
    ax.set(xticks=np.arange(indx1,indx2+1,(indx2-indx1)/2))
    ax.set(xticklabels=np.arange(indx1,indx2+1,(indx2-indx1)/2)*1000/fs)
    simpleaxis(ax)
    plt.ylim([yrange[0]*100,yrange[1]*100])
    ax.set(yticks=[])
    
    ax=plt.subplot(122)  
    r = abs(csd).max()/5
    im = ax.imshow(csd,origin='image',
                   vmin=-r, 
                   vmax=r, cmap='jet',
                   interpolation='nearest',aspect='auto')
    cb = plt.colorbar(im, ax=ax,orientation="vertical")
    #cb.set_label('CSD (A/m**3)')
    plt.xlabel('Time (ms)')  
    simpleaxis(ax)
    plt.xlim([indx1,indx2+1])
    ax.set(xticks=np.arange(indx1,indx2+1,(indx2-indx1)/2))
    ax.set(xticklabels=np.arange(indx1,indx2+1,(indx2-indx1)/2)*1000/fs)
    plt.ylim([yrange[0],yrange[1]])
    ax.set(yticks=[])
        
    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.4,
               wspace=0.4)
    return fig



    ax.spines['left'].set_visible(False)

