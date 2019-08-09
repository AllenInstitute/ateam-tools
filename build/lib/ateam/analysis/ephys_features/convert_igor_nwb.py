#!/usr/bin/python
import os
import h5py
import sys
import shutil
import traceback
import subprocess
import numpy as np
import re
import os

import nwb
#from allensdk.internal.core.lims_pipeline_module import PipelineModule

# development/debugging code
#infile = "Ndnf-IRES2-dgCre_Ai14-256189.05.01-compressed.nwb"
#outfile = "foo.nwb"
#if len(sys.argv) == 1:
#    sys.argv.append(infile)
#    sys.argv.append(outfile)

# this script is meant to clone the core functionality of the 
#   existing (Igor) Hdf5->Nwb converter.
# the previous converter performed two distinct tasks. In this iteration,
#   those tasks will be split into separate modules. This module will
#   perform the file conversion. A second module will analyze the file
#   and extract sweep data


# show the nwb file list
nwb_file_list=list()

for file in os.listdir("./"):
    if file.endswith(".nwb"):
        nwb_file_list.append(os.path.join(file))

print 'Here are the abf files'
print nwb_file_list
print '\n'




# window for leading test pulse, in seconds
PULSE_LEN = 0.05
EXPERIMENT_START_TIME = 0.1


def main(file_name):
    #module = PipelineModule()
    #jin = module.input_data()

# WRITE DOWN THE FILE NAMES HERE

    #infile = "H17.06.015.21.09.06.nwb"
    #outfile = "H17.06.015.21.09.06_converted.nwb"

    # read the convert the first nwb file
#    nwb_file_list=list()

#    for file in os.listdir("./"):
#        if file.endswith(".nwb"):
#            nwb_file_list.append(os.path.join(file))

    infile = file_name
    outfile = str(infile[0:-4]) +str('_converted') +str('.nwb')


    # a temporary nwb file must be created. this is that file's name
    tmpfile = outfile + ".tmp"

    # create temp file and make modifications to it using h5py
    shutil.copy2(infile, tmpfile)
    f = h5py.File(tmpfile, "a")
    # change dataset names in acquisition time series to match that
    #   of existing ephys NWB files
    # also rescale the contents of 'data' fields to match the scaling
    #   in original files
    acq = f["acquisition/timeseries"]

    sweep_nums = []
    for k, v in acq.iteritems():
        # parse out sweep number
        try:
#            print 'Here is the num'
            num_str=str(k)
            #print num_str
            num = int(re.findall("\d+", num_str)[0])
        except:
            print("Error - unexpected sweep name encountered in IGOR nwb file")
            print("Sweep called: '%s'" % k)
            print("Expecting 5-digit sweep number between chars 5 and 9")
#            sys.exit(1)
        swp = "Sweep_%d" % num
        # rename objects
        try:
            acq.move(k, swp)
            ts = acq[swp]
#           ts.move("stimulus_description", "aibs_stimulus_description")        


        except:
            print("*** Error renaming HDF5 object in %s" % swp)
            type_, value_, traceback_ = sys.exc_info()
            print traceback.print_tb(traceback_)
            sys.exit(1)
        # rescale contents of data so conversion is 1.0
        try:
            data = ts["data"]
            scale = 1
            #float(data.attrs["conversion"])
            data[...] = data.value * scale
            data.attrs["conversion"] = 1.0
        except:
            print("*** Error rescaling data in %s" % swp)
            type_, value_, traceback_ = sys.exc_info()
            print traceback.print_tb(traceback_)
            sys.exit(1)
        # keep track of sweep numbers
        sweep_nums.append("%d"%num)
        
    ###################################
    #... ditto for stimulus time series

    stim = f["stimulus/presentation"]

    for k, v in stim.iteritems():
        # parse out sweep number
        try:
#            print 'Here is the num'
            num_str=str(k)
            #print num_str
            num = int(re.findall("\d+", num_str)[0])

        except:
            print("Error - unexpected sweep name encountered in IGOR nwb file")
            print("Sweep called: '%s'" % k)
            print("Expecting 5-digit sweep number between chars 5 and 9")
#           sys.exit(1)
        swp = "Sweep_%d" % num
        try:
            stim.move(k, swp)
        except:
            print("Error renaming HDF5 group from %s to %s" % (k, swp))
#            sys.exit(1)
        # rescale contents of data so conversion is 1.0
    #try:
        ts = stim[swp]
        data = ts["data"]
        scale = 1
        #float(data.attrs["conversion"])
        data[...] = data.value * scale
        data.attrs["conversion"] = 1.0            
        # This part assigns name according to the sweep numbers
        stimulus_path=str('acquisition/timeseries/Sweep_')+str(num)+str('/aibs_stimulus_description/')

        print stimulus_path

        #stimulus_name= str(f.get(stimulus_path).value)
        stimulus_name = str(f.get(stimulus_path).value)
                
        # Rename the stimuli according to the names in the optimization
        

        # TEST STIMULI

        if 'EXTPGGAEND141203' in stimulus_name:
            stimulus_name ='Test stimuli'

        if 'EXTPEXPEND141203' in stimulus_name:
            stimulus_name ='Test stimuli'

        if 'EXTPBLWOUT141203' in stimulus_name:
            stimulus_name ='Test stimuli'

        if 'EXTPCllATT141203' in stimulus_name:
            stimulus_name ='Test stimuli'

        if 'EXTPRSCHEK150209' in stimulus_name:
            stimulus_name ='Test stimuli'

        if 'EXTPSMOKET141203' in stimulus_name:
            stimulus_name ='Test stimuli'

        if 'EXTPSAFETY141203' in stimulus_name:
            stimulus_name ='Test stimuli'

        if 'EXTPINBATH141203' in stimulus_name:
            stimulus_name ='Test stimuli'

        if 'EXTPBREAKN141203' in stimulus_name:
            stimulus_name ='Test stimuli'


        # SHORT SQUARE STIMULI

        if 'C2SSHM80CS150112' in stimulus_name:
            stimulus_name ='Short Square'

        if 'C2SSHM80FN150112' in stimulus_name:
            stimulus_name ='Short Square'

        if 'C2SSHM80FN141203' in stimulus_name:
            stimulus_name ='Short Square'

        if 'C2SSHM80CS141203' in stimulus_name:
            stimulus_name ='Short Square'


        # LONG DC STIMULI

        if 'C1LSCOARSE150216_DA_0' in stimulus_name:
            stimulus_name ='Long Square'

        if 'C1LSFINEST150112_DA_0' in stimulus_name:
            stimulus_name ='Long Square'

        if 'C1LSFINEST150112' in stimulus_name:
            stimulus_name ='Long Square'    

        if 'C1NSSEED_1150216' in stimulus_name:
            stimulus_name ='Long Square'                

        if 'C1LSCRESEP170313_DA_0' in stimulus_name:
            stimulus_name ='Long Square'

        if 'C1LSCOARSE150112' in stimulus_name:
            stimulus_name ='Long Square'

        if 'C1LSCOARSEMICRO' in stimulus_name:
            stimulus_name ='Long Square'

        if 'C1LSCOARSE150216' in stimulus_name:
            stimulus_name ='Long Square'

        if 'C1LSCOARSE150217' in stimulus_name:
            stimulus_name ='Long Square'

        if 'C1LSFINESTMICRO' in stimulus_name:
            stimulus_name ='Long Square'

        if 'C1LSCOARSE141203' in stimulus_name:
            stimulus_name ='Long Square'

        if 'C1LSFINEST150112' in stimulus_name:
            stimulus_name ='Long Square'

        if 'C1LSFINEST141203' in stimulus_name:
            stimulus_name ='Long Square'


		# RAMP stimuli

        if 'C1RP25PR1S141203' in stimulus_name:
            stimulus_name ='Ramp'

        # Noise stimuli

        if 'C1NSSEED_1150216' in stimulus_name:
            stimulus_name ='Noise 1'

        if 'C1NSSEED_1150112' in stimulus_name:
            stimulus_name ='Noise 1'

        if 'C1NSSEED_1141203' in stimulus_name:
            stimulus_name ='Noise 1'

        if 'C1NSSEED_1150217' in stimulus_name:
            stimulus_name ='Noise 1'


        if num<=4:
            stimulus_name ='Test stimuli'


        print stimulus_name

        
        # CHECK THE VOLTAGE
        voltage_path=str("acquisition/timeseries/Sweep_") + str(num) + str('/data')
        voltage = f.get(voltage_path).value*1e3 # get the mV
        # set the threshold value in mV

        # check for stim amplitude to get rid of strange sweeps
#        if max(abs(voltage*1e3))<1:
#            stimulus_name = 'Test stimuli'

        # special exception for sweeps
        #if max(voltage*1e3) < 0.001:
        #    if max(voltage*1e3) > -0.001:
        #        stimulus_name = 'Test stimuli'

        ts.create_dataset("aibs_stimulus_name",data=stimulus_name)

        # record the start of the stimulus in sec            
        ts.create_dataset('stimulus_start_sec',data=float(0+EXPERIMENT_START_TIME))

        # record the end of the stimulus in sec
        rate_path=str("acquisition/timeseries/Sweep_") + str(num) + str('/starting_time')
        sampling_rate = f.get(rate_path).attrs['rate']
        
        # show the sampling rate
#        print
#        print 'Sampling rate'
#        print sampling_rate

        samples_path=str("stimulus/presentation/Sweep_") + str(num) + str('/num_samples')
        n_samples=f.get(samples_path).value        
        ts.create_dataset('stimulus_end_sec',data=float(n_samples/sampling_rate))        

         
        #ts.create_dataset("aibs_stimulus_name",data=stimulus_name)

        # record the maximal stimulus ampitude
        stimulus_path=str("stimulus/presentation/Sweep_") + str(num) + str('/data')
        stimulus=f.get(stimulus_path).value
        # need to start searching for the maximal value after the test step!

        # get the stimulus amplitude based on gradient
        gradient_f = np.gradient(stimulus)
        signal_max = max(np.gradient(stimulus))
        signal_min = min(np.gradient(stimulus))

        # find the max/min of the gradient
        first_ind = np.where(gradient_f == signal_max)[0][0]
        second_ind = np.where(gradient_f == signal_min)[0][0]

        # check for the first and second indexes
        if first_ind > second_ind:
            start_ind = second_ind
            end_ind = first_ind
        elif first_ind < second_ind:
            start_ind = first_ind
            end_ind = second_ind


        if stimulus_name == 'Long Square':
            stimulus_amp = np.mean(stimulus[start_ind:end_ind])
#            print 'stimulus_amp'
#            print stimulus_amp

        if stimulus_name == 'Ramp':
            stimulus_amp = np.mean(stimulus[end_ind])
#            print 'stimulus_amp'
#            print stimulus_amp


        if stimulus_name == 'Test stimuli':
            stimulus_amp = 0

#        print 'Stimulus amplitude' + ' Sweep:' + str(num)
#        print float(stimulus_amp*1e12)
#        print

        #if np.min(stimulus)<0:
        #    stimulus_peak=np.min(stimulus)
        #else:
        #    stimulus_peak=np.max(stimulus)

        # stimulus amplitude in pA

        ts.create_dataset('stimulus_amplitude_pa',data=stimulus_amp*1e12)

        # sweep epoch
        #t0 = ts["starting_time"].value
        #rate = float(ts["starting_time"].attrs["rate"])

#except:
        #print("*** Error rescaling data in %s" % swp)
    #    type_, value_, traceback_ = sys.exc_info()
    #    print traceback.print_tb(traceback_)
    #    sys.exit(1)
        
    f.close()

    ####################################################################
    # re-open file w/ nwb library and add indexing (epochs)
    
    '''
    nd = nwb.NWB(filename=tmpfile, modify=True)
    for num in sweep_nums:
        ts = nd.file_pointer["acquisition/timeseries/Sweep_" + num]
        # sweep epoch
        t0 = ts["starting_time"].value
        rate = float(ts["starting_time"].attrs["rate"])
        n = float(ts["num_samples"].value)
        t1 = t0 + (n-1) * rate
        ep = nd.create_epoch("Sweep_" + num, t0, t1)
        
        ep.add_timeseries("stimulus", "stimulus/presentation/Sweep_"+num)
        ep.add_timeseries("response", "acquisition/timeseries/Sweep_"+num)
        ep.finalize()
        if "CurrentClampSeries" in ts.attrs["ancestry"]:
            # test pulse epoch
            t0 = ts["starting_time"].value
            t1 = t0 + PULSE_LEN
            ep = nd.create_epoch("TestPulse_" + num, t0, t1)
            ep.add_timeseries("stimulus", "stimulus/presentation/Sweep_"+num)
            ep.add_timeseries("response", "acquisition/timeseries/Sweep_"+num)
            ep.finalize()
            # experiment epoch
            t0 = ts["starting_time"].value
            t1 = t0 + (n-1) * rate
            t0 += EXPERIMENT_START_TIME
            ep = nd.create_epoch("Experiment_" + num, t0, t1)
            ep.add_timeseries("stimulus", "stimulus/presentation/Sweep_"+num)
            ep.add_timeseries("response", "acquisition/timeseries/Sweep_"+num)
            ep.finalize()
    nd.close()

    '''

    # rescaling the contents of the data arrays causes the file to grow
    # execute hdf5-repack to get it back to its original size
    try:
        print("Repacking hdf5 file with compression")
        process = subprocess.Popen(["h5repack", "-f", "GZIP=4", tmpfile, outfile], stdout=subprocess.PIPE)
        process.wait()
    except:
        print("Unable to run h5repack on temporary nwb file")
        print("--------------------------------------------")
        raise

    try:
        print("Removing temporary file")
        os.remove(tmpfile)
    except:
        print("Unable to delete temporary file ('%s')" % tmpfile)
        raise

    # done (nothing to return)
    #module.write_output_data({})
    return


# show the nwb file list
nwb_file_list=list()

for file in os.listdir("./"):
    if file.endswith(".nwb"):
        nwb_file_list.append(os.path.join(file))

print 'Here are the abf files'
print nwb_file_list
print '\n'


for i in range(len(nwb_file_list)):
    main(nwb_file_list[i])


#if __name__=='__main__': main()

