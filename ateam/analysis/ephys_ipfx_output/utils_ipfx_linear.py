import json
import os
import numpy as np
import matplotlib.pyplot as plt
#import burst as brst


def convert_keys_to_string(dictionary):
    """Recursively converts dictionary keys to strings."""
    if not isinstance(dictionary, dict):
        return dictionary
    return dict((str(k), convert_keys_to_string(v)) 
        for k, v in dictionary.items())



def parse_json_ipfx(cell_id):
    
    # this function loads the ipfx generated json file, extracts features
    # and returns the output in the form of dictionary
    
    # open file if it exists
    
    print cell_id

    try:

        with open(cell_id + '/' + 'pipeline_output.json') as f:
            my_data = json.load(f)

        cell_record_features=my_data['feature_extraction']['cell_record']


        # remove setup-specific ephys features
        del cell_record_features['input_resistance_mohm']
        del cell_record_features['vm_for_sag']
        del cell_record_features['rheobase_sweep_num']
        del cell_record_features['thumbnail_sweep_num']
        del cell_record_features['electrode_0_pa']
        del cell_record_features['input_access_resistance_ratio']
        del cell_record_features['blowout_mv']
        del cell_record_features['seal_gohm']
        del cell_record_features['initial_access_resistance_mohm']


#        cell_record_features['burst_index_n']=np.nan
#        cell_record_features['burst_index_t']=np.nan
#        cell_record_features['burst_ratio']=np.nan
#        cell_record_features['burst_index_n_max']=np.nan
#        cell_record_features['burst_index_t_max']=np.nan
#        cell_record_features['burst_ratio_max']=np.nan

        # remove spike time related features
#        del cell_record_features['trough_t_ramp']
#        del cell_record_features['peak_t_ramp']
#        del cell_record_features['fast_trough_t_ramp']
#        del cell_record_features['slow_trough_t_ramp']
#        del cell_record_features['threshold_t_ramp']

#        del cell_record_features['trough_t_long_square']
#        del cell_record_features['peak_t_long_square']
#        del cell_record_features['fast_trough_t_long_square']
#        del cell_record_features['slow_trough_t_long_square']
#        del cell_record_features['threshold_t_long_square']

#        del cell_record_features['trough_t_short_square']
#        del cell_record_features['peak_t_short_square']
#        del cell_record_features['fast_trough_t_short_square']
#        del cell_record_features['slow_trough_t_short_square']
#        del cell_record_features['threshold_t_short_square']

        # remove short square related features (there is an error here)
#        del cell_record_features['peak_v_short_square']
#        del cell_record_features['upstroke_downstroke_ratio_short_square']
#        del cell_record_features['threshold_v_short_square']
#        del cell_record_features['threshold_i_short_square']
#        del cell_record_features['fast_trough_v_short_square']
#        del cell_record_features['slow_trough_v_short_square']
#        del cell_record_features['trough_v_short_square']


# add additional bursting features

        '''
        try:
        	import burst as brst
        	burst_dict=brst.process_cell(cell_id, use_silence=True)

        	# mean characteristics
        	burst_index_n=0
        	burst_index_t=0
        	burst_ratio=0

        	# max characteristics
        	burst_index_n_max=0
        	burst_index_t_max=0
        	burst_ratio_max=0

        	for i in np.arange(len(burst_dict)):
        	    # compute the mean
        	    burst_index_n=burst_index_n + burst_dict[i]['burst_index_n']
        	    burst_index_t=burst_index_t + burst_dict[i]['burst_index_t']
        	    burst_ratio=burst_ratio + burst_dict[i]['burst_ratio']

        	    # compute the max values
        	    if burst_index_n > burst_index_n_max:
        		burst_index_n_max = burst_index_n
        	    if burst_index_t > burst_index_n_max:
        		burst_index_t_max = burst_index_t
        	    if burst_ratio > burst_ratio:
        		burst_ratio_max = burst_ratio

        	# get the average values
        	burst_index_n=burst_index_n/i
        	burst_index_t=burst_index_t/i
        	burst_ratio=burst_ratio/i

        	# add features to the dictionary
        	cell_record_features['burst_index_n']=burst_index_n
        	cell_record_features['burst_index_t']=burst_index_t
        	cell_record_features['burst_ratio']=burst_ratio
        	cell_record_features['burst_index_n_max']=burst_index_n_max
        	cell_record_features['burst_index_t_max']=burst_index_t_max
        	cell_record_features['burst_ratio_max']=burst_ratio_max

        except (IOError, ValueError, RuntimeError, TypeError, NameError, KeyError):
            print 'Bursting metric problems'

        '''

        # create the dictionary of all sweeps
        sweep_len=len(my_data['sweep_extraction']['sweep_features'])

        # features dict from json file
        features_sweeps_dict={}

        for sweeps in np.arange(sweep_len):

            # process the sweep number
            sweep_number=my_data['sweep_extraction']['sweep_features'][sweeps]['sweep_number']

            # process the sweep type
            stimulus_name=str(my_data['sweep_extraction']['sweep_features'][sweeps]['stimulus_name'])


            # check the stimulus type
            if stimulus_name == 'Long Square':

                current_sweep_dict=my_data['feature_extraction']['sweep_features'][str(sweep_number)]['spikes']

                # dictionary for all sweeps
                features_sweeps_dict[sweep_number]=current_sweep_dict

                features_sweeps_dict[sweep_number]


        # cell features dictionary for scaling features
        cell_features={}


        # lists for corresponding feature values
        downstroke=[]
        downstroke_index=[]
        downstroke_t=[]
        downstroke_v=[]
        fast_trough_index=[]
        fast_trough_t=[]
        fast_trough_v=[]
        peak_index=[]
        peak_t=[]
        peak_v=[]
        slow_trough_index=[]
        slow_trough_t=[]
        slow_trough_v=[]
        threshold_index=[]
        threshold_t=[]
        threshold_v=[]
        trough_index=[]
        trough_t=[]
        trough_v=[]
        upstroke=[]
        upstroke_downstroke_ratio=[]
        upstroke_index=[]
        upstroke_t=[]
        upstroke_v=[]
        width=[]


        # lists for corresponding current values
        downstroke_current=[]
        downstroke_index_current=[]
        downstroke_t_current=[]
        downstroke_v_current=[]
        fast_trough_index_current=[]
        fast_trough_t_current=[]
        fast_trough_v_current=[]
        peak_index_current=[]
        peak_t_current=[]
        peak_v_current=[]
        slow_trough_index_current=[]
        slow_trough_t_current=[]
        slow_trough_v_current=[]
        threshold_index_current=[]
        threshold_t_current=[]
        threshold_v_current=[]
        trough_index_current=[]
        trough_t_current=[]
        trough_v_current=[]
        upstroke_current=[]
        upstroke_downstroke_ratio_current=[]
        upstroke_index_current=[]
        upstroke_t_current=[]
        upstroke_v_current=[]
        width_current=[]


        # extract values from the dictionaries
        for sweep_name in features_sweeps_dict:

        # print the value of the first spike if there is one
            if features_sweeps_dict[sweep_name]:

                # record all features associated with the first spike

                # downstroke
                if features_sweeps_dict[sweep_name][0]['downstroke'] is not None:
                    downstroke.append(features_sweeps_dict[sweep_name][0]['downstroke'])
                    downstroke_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # downstroke_index
                if features_sweeps_dict[sweep_name][0]['downstroke_index'] is not None:
                    downstroke_index.append(features_sweeps_dict[sweep_name][0]['downstroke_index'])
                    downstroke_index_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # downstroke voltage time
                if features_sweeps_dict[sweep_name][0]['downstroke_t'] is not None:
                    downstroke_t.append(features_sweeps_dict[sweep_name][0]['downstroke_t'])
                    downstroke_t_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # downstroke voltage
                if features_sweeps_dict[sweep_name][0]['downstroke_v'] is not None:
                    downstroke_v.append(features_sweeps_dict[sweep_name][0]['downstroke_v'])
                    downstroke_v_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # fast_trough_index
                if features_sweeps_dict[sweep_name][0]['fast_trough_index'] is not None:
                    fast_trough_index.append(features_sweeps_dict[sweep_name][0]['fast_trough_index'])
                    fast_trough_index_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # fast_trough_index time
                if features_sweeps_dict[sweep_name][0]['fast_trough_t'] is not None:
                    fast_trough_t.append(features_sweeps_dict[sweep_name][0]['fast_trough_t'])
                    fast_trough_t_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # fast_trough_v
                if features_sweeps_dict[sweep_name][0]['fast_trough_v'] is not None:
                    fast_trough_v.append(features_sweeps_dict[sweep_name][0]['fast_trough_v'])
                    fast_trough_v_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # peak_index
                if features_sweeps_dict[sweep_name][0]['peak_index'] is not None:
                    peak_index.append(features_sweeps_dict[sweep_name][0]['peak_index'])
                    peak_index_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # peak_t
                if features_sweeps_dict[sweep_name][0]['peak_t'] is not None:
                    peak_t.append(features_sweeps_dict[sweep_name][0]['peak_t'])
                    peak_t_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # peak_v
                if features_sweeps_dict[sweep_name][0]['peak_v'] is not None:
                    peak_v.append(features_sweeps_dict[sweep_name][0]['peak_v'])
                    peak_v_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # slow_trough_index
                if features_sweeps_dict[sweep_name][0]['slow_trough_index'] is not None:
                    slow_trough_index.append(features_sweeps_dict[sweep_name][0]['slow_trough_index'])
                    slow_trough_index_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # slow_trough_t
                if features_sweeps_dict[sweep_name][0]['slow_trough_t'] is not None:
                    slow_trough_t.append(features_sweeps_dict[sweep_name][0]['slow_trough_t'])
                    slow_trough_t_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # slow_trough_v
                if features_sweeps_dict[sweep_name][0]['slow_trough_v'] is not None:
                    slow_trough_v.append(features_sweeps_dict[sweep_name][0]['slow_trough_v'])
                    slow_trough_v_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # threshold_index
                if features_sweeps_dict[sweep_name][0]['threshold_index'] is not None:
                    threshold_index.append(features_sweeps_dict[sweep_name][0]['threshold_index'])
                    threshold_index_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # threshold_t
                if features_sweeps_dict[sweep_name][0]['threshold_t'] is not None:
                    threshold_t.append(features_sweeps_dict[sweep_name][0]['threshold_t'])
                    threshold_t_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # threshold_v
                if features_sweeps_dict[sweep_name][0]['threshold_v'] is not None:
                    threshold_v.append(features_sweeps_dict[sweep_name][0]['threshold_v'])
                    threshold_v_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # trough_index
                if features_sweeps_dict[sweep_name][0]['trough_index'] is not None:
                    trough_index.append(features_sweeps_dict[sweep_name][0]['trough_index'])
                    trough_index_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # trough_t
                if features_sweeps_dict[sweep_name][0]['trough_t'] is not None:
                    trough_t.append(features_sweeps_dict[sweep_name][0]['trough_t'])
                    trough_t_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # trough_v
                if features_sweeps_dict[sweep_name][0]['trough_v'] is not None:
                    trough_v.append(features_sweeps_dict[sweep_name][0]['trough_v'])
                    trough_v_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # upstroke
                if features_sweeps_dict[sweep_name][0]['upstroke'] is not None:
                    upstroke.append(features_sweeps_dict[sweep_name][0]['upstroke'])
                    upstroke_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # upstroke_downstroke_ratio
                if features_sweeps_dict[sweep_name][0]['upstroke_downstroke_ratio'] is not None:
                    upstroke_downstroke_ratio.append(features_sweeps_dict[sweep_name][0]['upstroke_downstroke_ratio'])
                    upstroke_downstroke_ratio_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # upstroke_index
                if features_sweeps_dict[sweep_name][0]['upstroke_index'] is not None:
                    upstroke_index.append(features_sweeps_dict[sweep_name][0]['upstroke_index'])
                    upstroke_index_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # upstroke_t
                if features_sweeps_dict[sweep_name][0]['upstroke_t'] is not None:
                    upstroke_t.append(features_sweeps_dict[sweep_name][0]['upstroke_t'])
                    upstroke_t_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # upstroke_v
                if features_sweeps_dict[sweep_name][0]['upstroke_v'] is not None:
                    upstroke_v.append(features_sweeps_dict[sweep_name][0]['upstroke_v'])
                    upstroke_v_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])

                # width
                if features_sweeps_dict[sweep_name][0]['width'] is not None:
                    width.append(features_sweeps_dict[sweep_name][0]['width'])
                    width_current.append(features_sweeps_dict[sweep_name][0]['peak_i'])


        # downstroke
        if downstroke_current:
            x=downstroke_current
            y=downstroke
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['downstroke_scaling']=m
            cell_features['downstroke_rheobase']=y[x.index(min(x))]


        # downstroke_index
        if downstroke_index_current:
            x=downstroke_index_current
            y=downstroke_index
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['downstroke_index_scaling']=m
            cell_features['downstroke_index_rheobase']=y[x.index(min(x))]


        # downstroke_t
        if downstroke_t_current:
            x=downstroke_t_current
            y=downstroke_t
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['downstroke_t_scaling']=m
            cell_features['downstroke_t_rheobase']=y[x.index(min(x))]


        # downstroke_v
        if downstroke_v_current:
            x=downstroke_v_current
            y=downstroke_v
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['downstroke_v_scaling']=m
            cell_features['downstroke_v_rheobase']=y[x.index(min(x))]


        # fast_trough_index
        if fast_trough_index_current:
            x=fast_trough_index_current
            y=fast_trough_index
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['fast_trough_index_scaling']=m
            cell_features['fast_trough_index_rheobase']=y[x.index(min(x))]


        # fast_trough_t
        if fast_trough_t_current:
            x=fast_trough_t_current
            y=fast_trough_t
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['fast_trough_t_index_scaling']=m
            cell_features['fast_trough_t_rheobase']=y[x.index(min(x))]


        # fast_trough_v
        if fast_trough_v_current:
            x=fast_trough_v_current
            y=fast_trough_v
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['fast_trough_v_index_scaling']=m
            cell_features['fast_trough_v_rheobase']=y[x.index(min(x))]


        # peak_index
        if peak_index_current:
            x=peak_index_current
            y=peak_index
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['peak_index_index_scaling']=m
            cell_features['peak_index_rheobase']=y[x.index(min(x))]


        # peak_t
        if peak_t_current:
            x=peak_t_current
            y=peak_t
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['peak_t_index_scaling']=m
            cell_features['peak_t_rheobase']=y[x.index(min(x))]


        # peak_v
        if peak_v_current:
            x=peak_v_current
            y=peak_v
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['peak_v_index_scaling']=m
            cell_features['peak_v_rheobase']=y[x.index(min(x))]


        # slow_trough_index
        if slow_trough_index_current:
            x=slow_trough_index_current
            y=slow_trough_index
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['slow_trough_index_scaling']=m
            cell_features['slow_trough_index_rheobase']=y[x.index(min(x))]


        # slow_trough_t
        if slow_trough_t_current:
            x=slow_trough_t_current
            y=slow_trough_t
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['slow_trough_t_scaling']=m
            cell_features['slow_trough_t_rheobase']=y[x.index(min(x))]


        # slow_trough_v
        if slow_trough_v_current:
            x=slow_trough_v_current
            y=slow_trough_v
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['slow_trough_v_scaling']=m
            cell_features['slow_trough_v_rheobase']=y[x.index(min(x))]


        # threshold_index
        if threshold_index_current:
            x=threshold_index_current
            y=threshold_index
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['threshold_index_scaling']=m
            cell_features['threshold_index_rheobase']=y[x.index(min(x))]


        # threshold_t
        if threshold_t_current:
            x=threshold_t_current
            y=threshold_t
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['threshold_t_scaling']=m
            cell_features['threshold_t_rheobase']=y[x.index(min(x))]

        # threshold_v
        if threshold_v_current:
            x=threshold_v_current
            y=threshold_v
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['threshold_v_scaling']=m
            cell_features['threshold_v_rheobase']=y[x.index(min(x))]


        # trough_index
        if trough_index_current:
            x=trough_index_current
            y=trough_index
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['trough_index_scaling']=m
            cell_features['trough_index_rheobase']=y[x.index(min(x))]


        # trough_t
        if trough_t_current:
            x=trough_t_current
            y=trough_t
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['trough_t_scaling']=m
            cell_features['trough_t_rheobase']=y[x.index(min(x))]


        # trough_v
        if trough_v_current:
            x=trough_v_current
            y=trough_v
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['trough_v_scaling']=m
            cell_features['trough_v_rheobase']=y[x.index(min(x))]


        # upstroke
        if upstroke_current:
            x=upstroke_current
            y=upstroke
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['upstroke_scaling']=m
            cell_features['upstroke_rheobase']=y[x.index(min(x))]


        # upstroke_downstroke_ratio
        if upstroke_downstroke_ratio_current:
            x=upstroke_downstroke_ratio_current
            y=upstroke_downstroke_ratio
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['upstroke_downstroke_ratio_scaling']=m
            cell_features['upstroke_downstroke_ratio_rheobase']=y[x.index(min(x))]


        # upstroke_index
        if upstroke_index_current:
            x=upstroke_index_current
            y=upstroke_index
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['upstroke_index_scaling']=m
            cell_features['upstroke_index_rheobase']=y[x.index(min(x))]


        # upstroke_t
        if upstroke_t_current:
            x=upstroke_t_current
            y=upstroke_t
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['upstroke_t_scaling']=m
            cell_features['upstroke_t_rheobase']=y[x.index(min(x))]


        # upstroke_v
        if upstroke_v_current:
            x=upstroke_v_current
            y=upstroke_v
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['upstroke_v_scaling']=m
            cell_features['upstroke_v_rheobase']=y[x.index(min(x))]


        # width
        if width_current:
            x=width_current
            y=width
            A = np.vstack([x, np.ones(len(y))]).T
            m, c = np.linalg.lstsq(A, y)[0]
            # record slope and intercent
            cell_features['width_scaling']=m
            cell_features['width_rheobase']=y[x.index(min(x))]

        # merge two dictionaries
        z = cell_features.copy()

        # currently we use only scaling features!
        z.update(cell_record_features)

        # convert all keys to string
#       z=convert_keys_to_string(z)

        # return the final dictionary
        return z

    # if errors - return an empty dictionary
    except (IOError, ValueError, RuntimeError, TypeError, NameError, KeyError):
        z={}
        return z

