######################################################
# Authors: Fahimeh Baftizadeh
# Date created: 4/1/2018
######################################################

import argparse

from bmtk.analyzer.build_network_table import *
from bmtk.analyzer.network_utils import StimType
from bmtk.analyzer.network_utils import ModelType
from bmtk.analyzer.network_utils import InputType


"""
This script is a wrapper for reobase_analysis/build_table.py to allow running on command line.
Only input parsing/validation is done here.
One difference is that inputs may be given in print format instead of explicit format (10 instead of -0.01)
"""
# args
parser = argparse.ArgumentParser()
parser.add_argument("network_id", help="cell gid, e.g. 313862022, 314900022, or 320668879", type=int)
parser.add_argument("stim_type", help="stimulus type, e.g. 'dc', 'sin'")
parser.add_argument("model_type", help="biophysical model, e.g. 'perisomatic', 'all_active'")
parser.add_argument("input_type", help="extracellular and intracellular stim, e.g. 'extrastim', 'extrastim_intrastim'")

# optional
parser.add_argument("-i", "--inputs", help="input values. may provide many. may use normal or display form for values, e.g. -0.06 == 60. Default is 10,20,...,100", default=default_inputs,nargs='+', type=float)
parser.add_argument("-t", "--trial", help="trial number for output folders, default 0", type=int, default=0)
sargs = parser.parse_args()

inputs = sargs.inputs

# Deal with input format
if np.abs(int(inputs[0])) == inputs[0]:
    # assume in display form
    inputs = [-1 * i/1000000 for i in inputs]

# Validate stim type
if sargs.stim_type not in [e.value for e in StimType]:
    raise ValueError('Invalid stim type provided (be sure type has been added to StimType class)')

if sargs.input_type not in [e.value for e in InputType]:
    raise ValueError('Invalid input type provided (be sure type has been added to InputType class)')

# Validate model type
if sargs.model_type not in [e.value for e in ModelType]:
    raise ValueError('Invalid model type provided (be sure type has been added to ModelType class)')

# Run it!
#build_dc(sargs.cell_gid, inputs, sargs.input_type, sargs.stim_type, sargs.model_type, sargs.trial)
build_table(sargs.network_id, sargs.input_type, sargs.model_type, sargs.stim_type, inputs,  sargs.trial, include_spike_phase = True, saved_data= False, input_freq=None)
