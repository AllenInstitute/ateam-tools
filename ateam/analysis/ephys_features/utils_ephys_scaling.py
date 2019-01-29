# Allen Institute Software License - This software license is the 2-clause BSD
# license plus a third clause that prohibits redistribution for commercial
# purposes without further permission.
#
# Copyright 2019. Allen Institute. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# 3. Redistributions for commercial purposes are not permitted without the
# Allen Institute's written permission.
# For purposes of this license, commercial purposes is the incorporation of the
# Allen Institute's software into anything for which you will charge fees or
# other compensation. Contact terms@alleninstitute.org for commercial licensing
# opportunities.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.


# import standard libraries
import numpy as np
import os
import scipy
# import h5py
import h5py
# import feature extractor
import allensdk
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor



def analyse_nwb_files(path_to_files):
    
    nwb_file_list=list()

    for file in os.listdir(path_to_files):
        if file.endswith(".nwb"):
            nwb_file_list.append(os.path.join(file))

    # create a dictionary of dictionaries
    # to store computed features
    files_features={}
    files_features_AP={}

    # dictionary of dictionaries for firing-rate
    files_features_fI={}

    # dictionaries to record spikes
    voltage_spike=[]
    time_spike=[]

    # dictionaries to record f-I curves
    firing_frequency=[]
    input_current=[]


    for l in range(0,len(nwb_file_list)):

        # load nwb as h5py
        file_name=path_to_files+str('/')+str(nwb_file_list[l])
        f = h5py.File(file_name,'r')

        # get the list of time series, access by names
        list_timeseries=f.get('acquisition/timeseries/').keys()

        print 'Analyzed file'
        print nwb_file_list[l]
        print '\n'

        print 'There are'
        print len(list_timeseries)
        print 'Sweeps'

        ### Action potential properties: Volate deap after AP (fast through-index)
        

        # Extraction of features for all sweeps
        sweep_characteristics={}

        # load the file
        f = h5py.File(file_name,'r')

        for i in range(1,len(list_timeseries)):

            # prepare the feature drict
            features_dict={}

            voltage_name=str('acquisition/timeseries/Sweep_')+str(i)+str('/data')
            current_name=str('stimulus/presentation/Sweep_')+str(i)+str('/data')
            rate_name=str('acquisition/timeseries/Sweep_')+str(i)+str('/starting_time')
            stim_amp_name=str('stimulus/presentation/Sweep_')+str(i)+str('/stimulus_amplitude_pa')
            stim_name=str('stimulus/presentation/Sweep_')+str(i)+str('/aibs_stimulus_name')

            # in mV
            voltage=f[voltage_name].value*1e3    
            # in pA
            current=f[current_name].value*1e12        
            # extracts the sampling rate for a given recording in Hz
            # sampling rate
            sampling_rate=f[rate_name].attrs['rate']        
            
            # substact the baseline from the input current, the average of the last 50 ms
            last_50_ms_steps=int(0.05*sampling_rate)
            current=current - np.mean(current[-last_50_ms_steps:])
            
            # time is in sec
            time = np.arange(0, len(voltage)) * (1.0 / int(sampling_rate))

            # decimate downsampling step if there are spikes
#            if max(voltage) > 0:
            voltage=scipy.signal.decimate(voltage,5)
            current=scipy.signal.decimate(current,5)                    
            time=scipy.signal.decimate(time,5)

                    # check for normal spikes
            if max(voltage) < 100:
                if min(voltage) > - 150:
                    current_amp=f[stim_amp_name].value

                    # record voltage deflection if there is a negative amplitude
                    if current_amp<0:
                        gradient_current=np.gradient(current)
                        signal_max=max(np.gradient(current))
                        signal_min=min(np.gradient(current))
                        # find the first and second indexes of the current step
                        first_ind=np.where(gradient_current == signal_max)[0][0]
                        second_ind=np.where(gradient_current == signal_min)[0][0]
                        # check for the first and second indexes
                        if first_ind>second_ind:
                            start_ind=second_ind
                            end_ind=first_ind
                        elif first_ind<second_ind:
                            start_ind=first_ind
                            end_ind=second_ind
                        # get the average voltage in the middle of the sweep
                        voltage_negative_amp=np.mean(voltage[(int(first_ind+second_ind)/2-int(0.1*sampling_rate)):(int(first_ind+second_ind)/2+int(0.1*sampling_rate))])        

                    # calculates the voltage baseline at the end of each sweep, after application of the test pulse
                    voltage_baseline=np.mean(voltage[0:int(0.1*sampling_rate)])

                    # FITTING THE TIME CONSTANT
                    # the fit is based on current step injections only!
                    if f[stim_name].value == 'Long Square':
                        if max(voltage) <= 0:
                            gradient_current=np.gradient(current)
                            signal_max=max(np.gradient(current))
                            signal_min=min(np.gradient(current))
                            # find the first and second indexes of the current step
                            first_ind=np.where(gradient_current == signal_max)[0][0]
                            second_ind=np.where(gradient_current == signal_min)[0][0]
                            # check for the first and second indexes
                            if first_ind>second_ind:
                                start_ind=second_ind
                                end_ind=first_ind
                            elif first_ind<second_ind:
                                start_ind=first_ind
                                end_ind=second_ind

                            # check for zero stimulus amplitude, if not zero, then fit
                            stim_amp = abs(np.mean(current[start_ind:end_ind]))
                            if stim_amp > 10:            
                                voltage_to_fit=voltage[int(start_ind):int(np.round((start_ind + 0.05*sampling_rate)))]
                                # a, inv_tau, y0, we need inv_tau
                                taum=allensdk.ephys.ephys_features.fit_membrane_time_constant(voltage, time, time[int(start_ind)], time[int(np.round((start_ind + 0.05*sampling_rate)))])
                                    # add time const in ms to features
                                if taum[1]:
                                    features_dict['taum']=(1/taum[1])*1000

                    # calculate the timing of the input stimuli, start and end
                    gradient_current=np.gradient(current)
                    signal_max=max(np.gradient(current))
                    signal_min=min(np.gradient(current))
                    # find the first and second indexes of the current step
                    first_ind=np.where(gradient_current == signal_max)[0][0]
                    second_ind=np.where(gradient_current == signal_min)[0][0]
                    # check for the first and second indexes
                    if first_ind>second_ind:
                        start_ind=second_ind
                        end_ind=first_ind
                    elif first_ind<second_ind:
                        start_ind=first_ind
                        end_ind=second_ind
                    #record the start/end time to the dictionary
                    features_dict['stim_start']=time[start_ind]
                    features_dict['stim_end']=time[end_ind]


                    # EXTRACTING SPIKE FEATURES
                    if f[stim_name].value == 'Long Square':
                        sweep = EphysSweepFeatureExtractor(t=time, v=voltage, i=current, start=0, end=time[-1])
                        sweep.process_spikes()

                        # RECORD FEATURES INTO THE DICTIONARY OF DICTIONARIES
                        all_spike_features=sweep.spike_feature_keys()
                        for j in range(0,len(all_spike_features)):
                            features_dict[all_spike_features[j]]=sweep.spike_feature(all_spike_features[j])
                    
                    # record the current in pA into the dictionary
                    features_dict['current_amp']=current_amp
                    # record the voltage deflection in mV into the dictionary
                    
                    if current_amp<0:
                        features_dict['voltage_negative_deflection']=voltage_negative_amp        
                    # record voltage baseline
                    features_dict['voltage_base']=voltage_baseline

                    # Record features in a dict of dict for various sweeps
                    sweep_characteristics[str('Sweep_')+str(i)]=features_dict


                    if f[stim_name].value == 'Long Square':
                        # check if there are spikes at all
                        sweep_name=str('Sweep_')+str(i)
                        if 'threshold_index' in sweep_characteristics[sweep_name]:
                            if np.mean(sweep_characteristics[sweep_name]['peak_i'][0]) >0:   # only for positive currents                            
                                # save the sweep if there is one spike (two spike for test)
                                if len(sweep_characteristics[sweep_name]['threshold_index']) <= 2:                                
                                    # record the spike voltage and time trace per spike
                                    spike_ind_start=sweep_characteristics[sweep_name]['threshold_index'][0]-int(0.005*sampling_rate)
                                    spike_ind_end=sweep_characteristics[sweep_name]['threshold_index'][0]+int(0.005*sampling_rate)
                                    # minimal voltage - first ms after the start
                                    min_voltage=np.mean(voltage[spike_ind_start:(spike_ind_start+int(0.001*sampling_rate))])
                                    voltage_spike.append((voltage[spike_ind_start:spike_ind_end]-min_voltage))
    #                                min_time=min(time[spike_ind_start:spike_ind_end])
    #                                time_spike.append((time[spike_ind_start:spike_ind_end]-min_time)*1000)
                                    min_time=np.arange(0, len((voltage[spike_ind_start:spike_ind_end]-min_voltage))) * (1.0 / int(sampling_rate))
                                    time_spike.append(min_time*1000)



        # create the dictionary

        # store features of a file
        nwb_features={}

        ### Passive properties: membrane constant

        # Estimate the taum from single pulses

        taum_list=[]

        # record all tam estimations into the list
        for sweep_name in sweep_characteristics:
            if 'taum' in sweep_characteristics[sweep_name]:
                taum_list.append(sweep_characteristics[sweep_name]['taum'])

        taum_list=np.array(taum_list)

        # Show taum
        print 'Initial taum estimate'
        print taum_list
        print '\n'

        # get only taum between 5 and 100 ms
        taum_list=taum_list[np.where( taum_list < 100 )]
        taum_list=taum_list[np.where( taum_list > 5 )]

        print 'Outliers removed'
        print taum_list

        # compute the values
        mean = np.mean(taum_list)
        sigma = np.std(taum_list)

        print 'The average membrane time constant in ms'
        print np.round(mean,2)

        nwb_features['taum_mean']=mean
        nwb_features['taum_sigma']=sigma

        # plot the taum as a function of sweep number

        taum_list=[]
        sweep_number=[]

        # record all tam estimations into the list
        for sweep_name in sweep_characteristics:
                if 'taum' in sweep_characteristics[sweep_name]:
                # get only the small taum between 5 and 100 ms to remove artefacts        
                    if sweep_characteristics[sweep_name]['taum'] <100:            
                        if sweep_characteristics[sweep_name]['taum'] >5:
                            string = sweep_name
                            # just gets the first element from the one valued list
                            number=[int(s) for s in string.split('_') if s.isdigit()][0]

                            sweep_number.append(number)
                            # check for NaNs
                            if sweep_characteristics[sweep_name]['taum']:                                                        
                                taum_list.append(sweep_characteristics[sweep_name]['taum'])

        # plot the results

        ### Passive properties: average voltage

        #
        # plot the taum as a function of sweep number
        voltage_base_list=[]

        # record all tam estimations into the list
        for sweep_name in sweep_characteristics:

                # get only the small taum between 5 and 100 ms to remove artefacts        
        #        if sweep_characteristics[sweep_name]['taum'] <100:            
        #            if sweep_characteristics[sweep_name]['taum'] >5:
                        string = sweep_name                
                        voltage_base_list.append(sweep_characteristics[sweep_name]['voltage_base'])

        # compute the values
        mean = np.mean(voltage_base_list)
        sigma = np.std(voltage_base_list)

        print 'The average membrane potential in mV'
        print np.round(mean,2)

        nwb_features['voltage_base_mean']=mean
        nwb_features['voltage_base_sigma']=sigma


        ### Passive properties: input resistance

        # Calculate and plot the input resistance for negative sweeps

        negative_current=[]
        negative_voltage=[]

        for sweep_name in sweep_characteristics:
            # pick up only the responses with the negative voltage deflections
            if bool(sweep_characteristics[sweep_name].get('voltage_negative_deflection')) == True:
                negative_current.append(sweep_characteristics[sweep_name]['current_amp'])
                negative_voltage.append(sweep_characteristics[sweep_name]['voltage_negative_deflection'] -sweep_characteristics[sweep_name]['voltage_base'])


        # create arrays to fit data
        negative_current=np.array(negative_current)/1e3  # convert to nA
        negative_voltage=np.array(negative_voltage)

        # pick up only currents in a certain range to avoid negative fitting
        pos_idx=np.where((np.array(negative_current)>-0.300) & (np.array(negative_current)<0) & (np.array(negative_voltage)< -5 ) )
        negative_current=[negative_current[i] for i in pos_idx[0]]
        negative_voltage=[negative_voltage[i] for i in pos_idx[0]]
        
        # create a linear fit for the data
        if negative_current:
            if negative_voltage:
                A = np.vstack([negative_current, np.ones(len(negative_current))]).T
                m, c = np.linalg.lstsq(A, negative_voltage)[0]
                nwb_features['R_in']=m                
        else:
            nwb_features['R_in']=np.nan
            
        print 'Input resistance in MOhm'
        print nwb_features['R_in']

        ## Active properties

        ### Active properties: f-I curve

        # Plot the f-I curve

        # load the NWB file
        f = h5py.File(file_name,'r')

        #list of input currents
        stim_input=[]
        # list of resulting frequencies
        spike_freq=[]

        # Over the range of long square sweeps
        for sweep_name in sweep_characteristics:

            # check the stimulus for name
            stimulus_type_path=str('stimulus/presentation/') +str(sweep_name) +str('/aibs_stimulus_name')
            stimulus_type_name= str(f[stimulus_type_path].value)

            # INPUT PROCESSING
            if stimulus_type_name == 'Long Square':
                stimulus_path=str("stimulus/presentation/") +str(sweep_name) +str('/stimulus_amplitude_pa')
                stimulus_value=f[stimulus_path].value # to get the pA out of A
                # record the input current
                stim_input.append(stimulus_value)

                # SPIKES PROCESSING        
                if 'peak_t' in sweep_characteristics[sweep_name]:

                     stimulus_start_path=str("stimulus/presentation/") +str(sweep_name) +str('/starting_time')
                     all_spikes=sweep_characteristics[sweep_name]['peak_t']
                     num_spikes=len(all_spikes)
                     spike_freq.append(num_spikes/1) # all spikes during the trial

                # if there are no spikes, add 0
                else:
                    spike_freq.append(0)

        ### Active properties: f-I curve slope


        #data for the fit
        
        pos_idx=np.where((np.array(spike_freq)>0.5) & (np.array(stim_input)>0))
        
        input_value=[stim_input[i] for i in pos_idx[0]]
        feature_value=[spike_freq[i] for i in pos_idx[0]]

        #extract the theobase current and minimal frequency
        try:
            min_freq=min(feature_value)
            # get the value of the minimal element
            pos_min=np.argmin(feature_value)
            min_curr=input_value[pos_min]        
            # save the rheobase current and frequency
            nwb_features['Rheobase_current']=min_curr
            nwb_features['Rheobase_freq']=min_freq
        except (ValueError, RuntimeError, TypeError, NameError):
            nwb_features['Rheobase_current']=np.nan
            nwb_features['Rheobase_freq']=np.nan


        # REMOVE NON-UNIQUE VALUES FROM TWO LISTS!
        input_value_unique=[]
        feature_value_unique=[]    
        # convert them to arrays
        input_value_unique=np.array(input_value)
        feature_value_unique=np.array(feature_value)
        # rounding step
        input_value_unique=np.around(input_value_unique,decimals=0)

        # get only unique in input current elemnts
        input_value_unique, idx_unique = np.unique(input_value_unique, return_index=True)
        feature_value_unique=feature_value_unique[idx_unique]
        try:
            
            # correction for nA case
            if max(input_value_unique)<1:
                input_value_unique=input_value_unique*1e12

            # insert zero values until input current
            zero_steps=np.arange(0,min(input_value_unique)-1,1)
            input_value_unique=np.concatenate((zero_steps, input_value_unique), axis=0)
            zero_freq=np.zeros(len(zero_steps))
            feature_value_unique=np.concatenate((zero_freq, feature_value_unique), axis=0)

            # save the f-I curves data
            input_current.append(input_value_unique)
            firing_frequency.append(feature_value_unique)    


            #Make a linear fit
            A = np.vstack([input_value, np.ones(len(input_value))]).T

            m, c = np.linalg.lstsq(A, feature_value)[0]
            nwb_features['fI_slope']=m
            
        except (ValueError,RuntimeError, TypeError, NameError):
            nwb_features['fI_slope']=np.nan
                        

        # save fI slope


        ### Active properties: relative average peak voltage

        # Plot the average peak voltage as a function of the stimulus

        # load the NWB file
        f = h5py.File(file_name,'r')

        #list of instantaneous input currents
        input_value=[]
        # list of instantaneous voltages
        voltage_value=[]

        # Over the range of long square sweeps
        for sweep_name in sweep_characteristics:

            # check the stimulus for name
            stimulus_type_path=str('stimulus/presentation/') +str(sweep_name) +str('/aibs_stimulus_name')
            stimulus_type_name= str(f[stimulus_type_path].value)

            # currents processing
            if stimulus_type_name == 'Long Square':        
                # Record the threshold voltage for the first spike

                # check if there are spikes at all
                if 'peak_v' in sweep_characteristics[sweep_name]:
                    if np.mean(sweep_characteristics[sweep_name]['peak_i'][0]) >0:   # only for positive currents
                        if np.mean(sweep_characteristics[sweep_name]['peak_i'][0]) < 800:   # only for positive currents
                            try:
                                voltage_value.append(sweep_characteristics[sweep_name]['peak_v'][0]-sweep_characteristics[sweep_name]['threshold_v'][0])
                                input_value.append(sweep_characteristics[sweep_name]['peak_i'][0])
                            except:
                                print voltage_value
                                print input_value


        #Make a linear fit
        A = np.vstack([input_value, np.ones(len(input_value))]).T
        try:
            m, c = np.linalg.lstsq(A, voltage_value)[0]
            nwb_features['AP_height_slope']=m
            nwb_features['AP_height_rheobase']=voltage_value[input_value.index(min(input_value))]
        except (ValueError,RuntimeError, TypeError, NameError):
            nwb_features['AP_height_slope']=np.nan
            nwb_features['AP_height_rheobase']=np.nan




        print 'Input current in pA'
        print np.round(input_value,2)
        print '\n'
        print 'Average peak voltage'
        print np.round(voltage_value,2)
        print '\n'
        print 'Average peak voltage'
        print np.round(np.mean(voltage_value),2)

        # load the NWB file
        f = h5py.File(file_name,'r')

        #list of instantaneous input currents
        input_value=[]
        # list of instantaneous voltages
        width_value=[]

        # Over the range of long square sweeps
        for sweep_name in sweep_characteristics:

            # check the stimulus for name
            stimulus_type_path=str('stimulus/presentation/') +str(sweep_name) +str('/aibs_stimulus_name')
            stimulus_type_name= str(f[stimulus_type_path].value)

            # currents processing
            if stimulus_type_name == 'Long Square':        
                # Record the threshold voltage for the first spike
                if 'peak_v' in sweep_characteristics[sweep_name]:
                    if np.mean(sweep_characteristics[sweep_name]['peak_i'][0]) >0:   # only for positive currents
                        if np.mean(sweep_characteristics[sweep_name]['peak_i'][0]) < 800:   # only for positive currents
                            try:
                                if sweep_characteristics[sweep_name]['width'][0]:
                                    input_value.append(sweep_characteristics[sweep_name]['peak_i'][0])
                                    width_value.append(sweep_characteristics[sweep_name]['width'][0]*1000)
                            except (ValueError, RuntimeError, TypeError, NameError, IndexError):
                                print input_value
                                print width_value


        #Make a linear fit
        try:
            A = np.vstack([input_value, np.ones(len(width_value))]).T
            m, c = np.linalg.lstsq(A, width_value)[0]
            nwb_features['AP_width_slope']=m
            nwb_features['AP_width_rheobase']=width_value[input_value.index(min(input_value))]
        except (ValueError,RuntimeError, TypeError, NameError):
            nwb_features['AP_width_slope']=np.nan
            nwb_features['AP_width_rheobase']=np.nan




        print 'Input current in pA'
        print np.round(input_value,2)
        print '\n'
        print 'Average peak width in ms'
        print np.round(width_value,2)
        print '\n'
        print 'Average AP widtj trials in mV'
        print np.round(np.mean(width_value),2)
        print '\n'
        print 'Average threshold std in mV'
        print np.round(np.std(width_value),2)


        # export to eps
        #plt.savefig('W1_AP_width.eps', format='eps', dpi=300)



        ### Active properties: AP voltage treshold

        # Plot the average peak voltage as a function of the stimulus

        # load the NWB file
        f = h5py.File(file_name,'r')

        #list of instantaneous input currents
        input_value=[]
        # list of instantaneous voltages
        thr_value=[]

        # Over the range of long square sweeps
        for sweep_name in sweep_characteristics:

            # check the stimulus for name
            stimulus_type_path=str('stimulus/presentation/') +str(sweep_name) +str('/aibs_stimulus_name')
            stimulus_type_name= str(f[stimulus_type_path].value)

            # currents processing
            if stimulus_type_name == 'Long Square':        
                # Record the threshold voltage for the first spike
                if 'peak_v' in sweep_characteristics[sweep_name]:
                    if np.mean(sweep_characteristics[sweep_name]['peak_i'][0]) >0:   # only for positive currents
                        if np.mean(sweep_characteristics[sweep_name]['peak_i'][0]) < 800:   # only for positive currents
                            try:
                                thr_value.append(sweep_characteristics[sweep_name]['threshold_v'][0])
                                input_value.append(sweep_characteristics[sweep_name]['peak_i'][0])
                            except:
                                print thr_value
                                print input_value

        #Make a linear fit
        A = np.vstack([input_value, np.ones(len(thr_value))]).T
        try:
            m, c = np.linalg.lstsq(A, thr_value)[0]
            nwb_features['AP_thr_slope']=m
            nwb_features['AP_thr_rheobase']=thr_value[input_value.index(min(input_value))]
        except (ValueError,RuntimeError, TypeError, NameError):
            nwb_features['AP_thr_slope']=np.nan
            nwb_features['AP_thr_rheobase']=np.nan


        print 'Input current in pA'
        print np.round(input_value,2)
        print '\n'
        print 'AP threshold value in mV'
        print np.round(thr_value,2)
        print '\n'

        print 'Average start in mV'
        print np.round(np.mean(thr_value),2)
        print '\n'
        print 'Average start std in mV'
        print np.round(np.std(thr_value),2)

        # Plot the average peak voltage as a function of the stimulus

        # load the NWB file
        f = h5py.File(file_name,'r')

        #list of instantaneous input currents
        input_value=[]
        # list of instantaneous voltages
        feature_value=[]

        # Over the range of long square sweeps
        for sweep_name in sweep_characteristics:

            # check the stimulus for name
            stimulus_type_path=str('stimulus/presentation/') +str(sweep_name) +str('/aibs_stimulus_name')
            stimulus_type_name= str(f[stimulus_type_path].value)

            # currents processing
            if stimulus_type_name == 'Long Square':        
                # Record the threshold voltage for the first spike
                if 'peak_v' in sweep_characteristics[sweep_name]:
                    if np.mean(sweep_characteristics[sweep_name]['peak_i'][0]) >0:   # only for positive currents
                        if np.mean(sweep_characteristics[sweep_name]['peak_i'][0]) < 800:   # only for positive currents
                            try:
                                feature_value.append(sweep_characteristics[sweep_name]['fast_trough_v'][0])
                                input_value.append(sweep_characteristics[sweep_name]['peak_i'][0])
                            except (ValueError, RuntimeError, TypeError, NameError, IndexError):
                                print input_value
                                print feature_value

        #Make a linear fit
        A = np.vstack([input_value, np.ones(len(feature_value))]).T
        try:
            m, c = np.linalg.lstsq(A, feature_value)[0]
            nwb_features['AP_through_slope']=m
            nwb_features['AP_through_rheobase']=feature_value[input_value.index(min(input_value))]
        except (ValueError,RuntimeError, TypeError, NameError):
            nwb_features['AP_through_slope']=np.nan
            nwb_features['AP_through_rheobase']=np.nan


        print 'Input current in pA'
        print np.round(input_value,2)
        print '\n'
        print 'Voltage trough after AP'
        print np.round(feature_value,2)
        print '\n'

        print 'Mean in mV'
        print np.round(np.mean(feature_value),2)
        print '\n'
        print 'Std in mV'
        print np.round(np.std(feature_value),2)




        ### Action potential properties: AP voltage upstroke

        # Plot the average peak voltage as a function of the stimulus

        # load the NWB file
        f = h5py.File(file_name,'r')

        #list of instantaneous input currents
        input_value=[]
        # list of instantaneous voltages
        feature_value=[]

        # Over the range of long square sweeps
        for sweep_name in sweep_characteristics:

            # check the stimulus for name
            stimulus_type_path=str('stimulus/presentation/') +str(sweep_name) +str('/aibs_stimulus_name')
            stimulus_type_name= str(f[stimulus_type_path].value)

            # currents processing
            if stimulus_type_name == 'Long Square':        
                # Record the threshold voltage for the first spike
                if 'peak_v' in sweep_characteristics[sweep_name]:
                    if np.mean(sweep_characteristics[sweep_name]['peak_i'][0]) >0:   # only for positive currents
                        if np.mean(sweep_characteristics[sweep_name]['peak_i'][0]) < 800:   # only for positive currents
                            input_value.append(sweep_characteristics[sweep_name]['peak_i'][0])
                            feature_value.append(sweep_characteristics[sweep_name]['upstroke_v'][0])

        #Make a linear fit
        A = np.vstack([input_value, np.ones(len(feature_value))]).T
        try:
            m, c = np.linalg.lstsq(A, feature_value)[0]
            nwb_features['AP_upstroke_slope']=m
            nwb_features['AP_upstroke_rheobase']=feature_value[input_value.index(min(input_value))]
        except (ValueError,RuntimeError, TypeError, NameError):
            nwb_features['AP_upstroke_slope']=np.nan
            nwb_features['AP_upstroke_rheobase']=np.nan

            
        print 'Input current in pA'
        print np.round(input_value,2)
        print '\n'
        print 'Voltage trough after AP in mV/ms'
        print np.round(feature_value,2)
        print '\n'

        print 'Mean in mV/ms'
        print np.round(np.mean(feature_value),2)
        print '\n'
        print 'Std in mV/ms'
        print np.round(np.std(feature_value),2)


        ### Action potential properties: AP voltage downstroke

        # Plot the average peak voltage as a function of the stimulus

        # load the NWB file
        f = h5py.File(file_name,'r')

        #list of instantaneous input currents
        input_value=[]
        # list of instantaneous voltages
        feature_value=[]

        # Over the range of long square sweeps
        for sweep_name in sweep_characteristics:

            # check the stimulus for name
            stimulus_type_path=str('stimulus/presentation/') +str(sweep_name) +str('/aibs_stimulus_name')
            stimulus_type_name= str(f[stimulus_type_path].value)

            # currents processing
            if stimulus_type_name == 'Long Square':        
                # Record the threshold voltage for the first spike
                if 'peak_v' in sweep_characteristics[sweep_name]:
                    if np.mean(sweep_characteristics[sweep_name]['peak_i'][0]) >0:   # only for positive currents
                        if np.mean(sweep_characteristics[sweep_name]['peak_i'][0]) < 800:   # only for positive currents
                                try:
                                    feature_value.append(sweep_characteristics[sweep_name]['downstroke_v'][0])
                                    input_value.append(sweep_characteristics[sweep_name]['peak_i'][0])
                                except (ValueError, RuntimeError, TypeError, NameError, IndexError, KeyError):
                                    print feature_value
                                    print input_value


        #Make a linear fit
        A = np.vstack([input_value, np.ones(len(feature_value))]).T
        try:
            m, c = np.linalg.lstsq(A, feature_value)[0]
            nwb_features['AP_downstroke_slope']=m
            nwb_features['AP_downstroke_rheobase']=feature_value[input_value.index(min(input_value))]
        except (ValueError,RuntimeError, TypeError, NameError):
            nwb_features['AP_downstroke_slope']=np.nan
            nwb_features['AP_downstroke_rheobase']=np.nan


        print 'Input current in pA'
        print np.round(input_value,2)
        print '\n'
        print 'Voltage trough after AP in mV/ms'
        print np.round(feature_value,2)
        print '\n'

        print 'Mean in mV/ms'
        print np.round(np.mean(feature_value),2)
        print '\n'
        print 'Std in mV/ms'
        print np.round(np.std(feature_value),2)


        ### Action potential properties: AP volage upstroke/downstroke ratio

        # Plot the average peak voltage as a function of the stimulus

        # load the NWB file
        f = h5py.File(file_name,'r')

        #list of instantaneous input currents
        input_value=[]
        # list of instantaneous voltages
        feature_value=[]

        # Over the range of long square sweeps
        for sweep_name in sweep_characteristics:

            # check the stimulus for name
            stimulus_type_path=str('stimulus/presentation/') +str(sweep_name) +str('/aibs_stimulus_name')
            stimulus_type_name= str(f[stimulus_type_path].value)

            # currents processing
            if stimulus_type_name == 'Long Square':        
                # Record the threshold voltage for the first spike
                if 'peak_v' in sweep_characteristics[sweep_name]:
                    if np.mean(sweep_characteristics[sweep_name]['peak_i'][0]) >0:   # only for positive currents
                        if np.mean(sweep_characteristics[sweep_name]['peak_i'][0]) < 800:   # only for positive currents
                                try:
                                    feature_value.append(sweep_characteristics[sweep_name]['upstroke_downstroke_ratio'][0])
                                    input_value.append(sweep_characteristics[sweep_name]['peak_i'][0])
                                except (ValueError, RuntimeError, TypeError, NameError, IndexError, KeyError):
                                    print feature_value
                                    print input_value
                                

        #Make a linear fit
        try:
            A = np.vstack([input_value, np.ones(len(feature_value))]).T
            m, c = np.linalg.lstsq(A, feature_value)[0]
            nwb_features['AP_up/downstroke_slope']=m
            nwb_features['AP_up/downstroke_rheobase']=feature_value[input_value.index(min(input_value))]
        except (ValueError,RuntimeError, TypeError, NameError):
            nwb_features['AP_up/downstroke_slope']=np.nan
            nwb_features['AP_up/downstroke_rheobase']=np.nan

        print 'Input current in pA'
        print np.round(input_value,2)
        print '\n'
        print 'Voltage upstroke/downstroke ratio'
        print np.round(feature_value,2)
        print '\n'

        print 'Mean in []'
        print np.round(np.mean(feature_value),2)
        print '\n'
        print 'Std in []'
        print np.round(np.std(feature_value),2)

        ### Active properties, spike patterns: adaptation index

        # Plot the average adaptation index as a function of the input current

        # load the NWB file
        f = h5py.File(file_name,'r')

        #list of instantaneous input currents
        input_value=[]
        # list of instantaneous voltages
        index_value=[]

        # Over the range of long square sweeps
        for sweep_name in sweep_characteristics:

            # check the stimulus for name
            stimulus_type_path=str('stimulus/presentation/') +str(sweep_name) +str('/aibs_stimulus_name')
            stimulus_type_name= str(f[stimulus_type_path].value)

            # currents processing
            if stimulus_type_name == 'Long Square':        
                # Record the threshold voltage for the first spike
                if 'peak_v' in sweep_characteristics[sweep_name]:
                    if np.mean(sweep_characteristics[sweep_name]['peak_i'][0]) >0:   # only for positive currents
                        if np.mean(sweep_characteristics[sweep_name]['peak_i'][0]) < 800:   # only for positive currents
                                input_value.append(sweep_characteristics[sweep_name]['peak_i'][0])
                                isis=sweep_characteristics[sweep_name]['peak_t']          # grab AP peak time
                                isis=np.diff(isis)                                        # compute the inter-spike intervals
                                isis=allensdk.ephys.ephys_features.adaptation_index(isis) # compute the adaptation index
                                index_value.append(isis)                                  # record the index
                                
        # convert not index
        index_value=np.array(index_value)
        #get the nan index value
        nan_index=np.argwhere(np.isnan(index_value))
        # remove elements with the wrong index value
        index_value = index_value[np.isfinite(index_value)]
        # remove the corresponding element in input value if it is needed
        if len(nan_index)>0:
            input_value = [i for j, i in enumerate(input_value) if j not in nan_index]
            
        

        #Make a linear fit        
        x=np.array(input_value)
        y=np.array(index_value)

        # take only positive values
        x=x[np.where(y>0)]
        y=y[np.where(y>0)]


        #Make a linear fit
        A = np.vstack([x, np.ones(len(x))]).T
        try:
            m, c = np.linalg.lstsq(A, y)[0]
            # save the adaptation index
            nwb_features['adaptation_slope']=m
            nwb_features['adaptation_slope_rheobase']=y[np.argmin(x)]
        except (ValueError,RuntimeError, TypeError, NameError):
                nwb_features['adaptation_slope']=np.nan
                nwb_features['adaptation_slope_rheobase']=np.nan

                
        print 'Input current in pA'
        print np.round(input_value,2)
        print '\n'
        print 'Adaptation index value'
        print np.round(index_value,2)
        print '\n'
        print 'Average adaptation index over trials'
        print np.round(np.mean(index_value),2)


        ### Active properties, spike patterns: average interspike interval

        # Plot the average peak voltage as a function of the stimulus

        # load the NWB file
        f = h5py.File(file_name,'r')

        #list of instantaneous input currents
        input_value=[]
        # list of instantaneous voltages
        feature_value=[]

        # Over the range of long square sweeps
        for sweep_name in sweep_characteristics:

            # check the stimulus for name
            stimulus_type_path=str('stimulus/presentation/') +str(sweep_name) +str('/aibs_stimulus_name')
            stimulus_type_name= str(f[stimulus_type_path].value)

            # currents processing
            if stimulus_type_name == 'Long Square':        
                # Record the threshold voltage for the first spike
                if 'peak_v' in sweep_characteristics[sweep_name]:
                    if np.mean(sweep_characteristics[sweep_name]['peak_i'][0]) >0:   # only for positive currents
                        if np.mean(sweep_characteristics[sweep_name]['peak_i'][0]) < 800:   # only for positive currents
                                isis=sweep_characteristics[sweep_name]['peak_t']*1000          # grab AP peak time
                                isis=np.diff(isis)                                        # compute the inter-spike intervals
                                if len(isis) >=1:
                                    input_value.append(np.mean(sweep_characteristics[sweep_name]['peak_i'][0]))                
                                    feature_value.append(np.mean(isis))                       # compute the average value of ISIs

        #Make an exponential fit
        x=np.array(input_value)
        y=np.array(feature_value)
        
        try:
            coeff=np.polyfit(x, np.log(y), 1, w=np.sqrt(y))
            nwb_features['mean_ISI_slope']=coeff[0]
            nwb_features['mean_ISI_slope_rheobase']=feature_value[input_value.index(min(input_value))]
        except (ValueError,RuntimeError, TypeError, NameError):
            nwb_features['mean_ISI_slope']=np.nan
            nwb_features['mean_ISI_slope_rheobase']=np.nan


        print 'Input current in pA'
        print np.round(input_value,2)
        print '\n'
        print 'Average interspike interval in ms'
        print np.round(feature_value,2)
        print '\n'

        print 'Mean in []'
        print np.round(np.mean(feature_value),2)
        print '\n'
        print 'Std in []'
        print np.round(np.std(feature_value),2)


        ### Active properties, spike patterns: first ISI

        # Plot the average peak voltage as a function of the stimulus

        # load the NWB file
        f = h5py.File(file_name,'r')

        #list of instantaneous input currents
        input_value=[]
        # list of instantaneous voltages
        feature_value=[]

        # Over the range of long square sweeps
        for sweep_name in sweep_characteristics:

            # check the stimulus for name
            stimulus_type_path=str('stimulus/presentation/') +str(sweep_name) +str('/aibs_stimulus_name')
            stimulus_type_name= str(f[stimulus_type_path].value)

            # currents processing
            if stimulus_type_name == 'Long Square':        
                # Record the threshold voltage for the first spike
                if 'peak_v' in sweep_characteristics[sweep_name]:
                    if np.mean(sweep_characteristics[sweep_name]['peak_i'][0]) >0:   # only for positive currents
                        if np.mean(sweep_characteristics[sweep_name]['peak_i'][0]) < 800:   # only for positive currents
                                isis=sweep_characteristics[sweep_name]['peak_t']*1000          # grab AP peak time
                                isis=np.diff(isis)                                        # compute the inter-spike intervals
                                if len(isis) >=1:
                                    input_value.append(sweep_characteristics[sweep_name]['peak_i'][0])                
                                    feature_value.append(isis[0])                       # compute the average value of ISIs

        #Make an exponential fit
        x=np.array(input_value)
        y=np.array(feature_value)
        try:
            coeff=np.polyfit(x, np.log(y), 1, w=np.sqrt(y))
            nwb_features['First_ISI_slope']=coeff[0]
            nwb_features['First_ISI_slope_rheobase']=feature_value[input_value.index(min(input_value))]
        except (ValueError,RuntimeError, TypeError, NameError):
            nwb_features['First_ISI_slope']=np.nan
            nwb_features['First_ISI_slope_rheobase']=np.nan


        print 'Input current in pA'
        print np.round(input_value,2)
        print '\n'
        print 'Average interspike interval in ms'
        print np.round(feature_value,2)
        print '\n'

        print 'Mean in []'
        print np.round(np.mean(feature_value),2)
        print '\n'
        print 'Std in []'
        print np.round(np.std(feature_value),2)


        ### Active properties, spike patterns: first spike delay

        # Plot the average peak voltage as a function of the stimulus

        # load the NWB file
        f = h5py.File(file_name,'r')

        #list of instantaneous input currents
        input_value=[]
        # list of instantaneous voltages
        feature_value=[]

        # Over the range of long square sweeps
        for sweep_name in sweep_characteristics:

            # check the stimulus for name
            stimulus_type_path=str('stimulus/presentation/') +str(sweep_name) +str('/aibs_stimulus_name')
            stimulus_type_name= str(f[stimulus_type_path].value)

            # currents processing
            if stimulus_type_name == 'Long Square':        
                # Record the threshold voltage for the first spike
                if 'peak_v' in sweep_characteristics[sweep_name]:
                    if np.mean(sweep_characteristics[sweep_name]['peak_i'][0]) >0:   # only for positive currents
                        if np.mean(sweep_characteristics[sweep_name]['peak_i'][0]) < 800:   # only for positive currents
                            try:
                                feature_value.append((sweep_characteristics[sweep_name]['peak_t'][0]-sweep_characteristics[sweep_name]['stim_start'])*1000)               
                                isis=sweep_characteristics[sweep_name]['peak_t']*1000          # grab AP peak time
                                isis=np.diff(isis)                                        # compute the inter-spike intervals
                                input_value.append(sweep_characteristics[sweep_name]['peak_i'][0])
                            except:
                                print feature_value
                                print input_value
                      # compute the average value of ISIs
 

        # Make a linear fit
        x=np.array(input_value)
        y=np.array(feature_value)
        try:
            coeff=np.polyfit(x, np.log(y), 1, w=np.sqrt(y))
            nwb_features['Time_to_spike_slope']=coeff[0]    
            nwb_features['Time_to_spike_slope_rheobase']=feature_value[input_value.index(min(input_value))]
        except (ValueError,RuntimeError, TypeError, NameError):
            nwb_features['Time_to_spike_slope']=np.nan
            nwb_features['Time_to_spike_slope_rheobase']=np.nan


        print 'Input current in pA'
        print np.round(input_value,2)
        print '\n'
        print 'Time to 1st spike in ms'
        print np.round(feature_value,2)
        print '\n'

        print 'Mean in ms'
        print np.round(np.mean(feature_value),2)
        print '\n'
        print 'Std in ms'
        print np.round(np.std(feature_value),2)

        # add first spikes to the dictionary of dictionaries
        files_features_AP['first_spike_voltage']=voltage_spike
        files_features_AP['first_spike_time']=time_spike

        files_features_fI['input_current']=input_current
        files_features_fI['firing_frequency']=firing_frequency


        # record features into the dictionary of dictionaries
        files_features[nwb_file_list[l]]=nwb_features

        
    return files_features, files_features_AP, files_features_fI



# this function creates the list from files features

def matrix_from_dict(n_cells,n_features,label,files_features_dict):
    
    # +1 to add labels
    all_features_matrix=np.zeros((len(files_features_dict.keys()),n_features+1))

    # file_names list
    file_names=[]

    # counter for file names
    k=0
    
    for key in enumerate(files_features_dict.keys()):
        # save the filename
        file_names.append(key[1])
        
        # counter for features
        m=0
        
        # feature name list
        feature_names=[]
        
        for keys in enumerate(files_features_dict[key[1]]):
            # save features to the matrix
            all_features_matrix[k,0]=label
            m=m+1
            all_features_matrix[k,m]=files_features_dict[key[1]][keys[1]]
            feature_names.append(keys[1])
            
            # test code to debug the function
#            print 'key[1]'
#            print key[1]
#            print
#            print 'keys[1]'
#            print keys[1]

        k=k+1

    return all_features_matrix, file_names, feature_names
