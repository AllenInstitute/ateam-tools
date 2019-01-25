import numpy as np
import pandas as pd
import json

class BaseWaveform(object):
    """
    Abstraction of waveform class to ensure calculate method is implemented
    """

    def calculate(self, simulation_time):
        raise NotImplementedError("Implement specific waveform calculation")

class BaseWaveformType(object):
    """
    Specific waveform type
    """

    def __init__(self, waveform_config):
        self.amp      = float(waveform_config["amp"]) # units? mA?
        self.delay    = float(waveform_config["del"]) # ms
        self.duration = float(waveform_config["dur"]) # ms

    def is_active(self, simulation_time):
        stop_time = self.delay + self.duration
        return self.delay < simulation_time < stop_time



class WaveformTypeDC(BaseWaveformType, BaseWaveform):
    """
    DC (step) waveform
    """

    def __init__(self, waveform_config):
        super(WaveformTypeDC, self).__init__(waveform_config)

    def calculate(self, t): # TODO better name
        if self.is_active(t):
            return self.amp
        else:
            return 0

class WaveformTypeSin(BaseWaveformType, BaseWaveform):
    """
    Sinusoidal waveform
    """

    def __init__(self, waveform_config):
        super(WaveformTypeSin, self).__init__(waveform_config)
        self.freq         = float(waveform_config["freq"])   # Hz
        self.phase_offset = float(waveform_config.get("phase", np.pi))  # radians, optional
        self.amp_offset   = float(waveform_config.get("offset", 0)) # units? mA? optional

    def calculate(self, t): # TODO better name

        if self.is_active(t):
            f = self.freq / 1000. # Hz to mHz
            a = self.amp
            return a * np.sin(2 * np.pi * f * t + self.phase_offset) + self.amp_offset
        else:
            return 0


class WaveformCustom(BaseWaveform):
    """
    Custom waveform defined by csv file
    """

    def __init__(self, waveform_file, wf_dir):
        self.definition = pd.read_csv(wf_dir + waveform_file, sep='\t')

    def calculate(self, t):
        return np.interp(t, self.definition["time"], self.definition["amplitude"])


class ComplexWaveform(BaseWaveform):
    """
    Superposition of simple waveforms
    """
    def __init__(self, el_collection):
        self.electrodes = el_collection

    def calculate(self, t):
        val = 0
        for el in self.electrodes:
            val += el.calculate(t)

        return val


# mapping from 'shape' code to subclass, always lowercase
shape_classes = {
    'dc': WaveformTypeDC,
    'sin': WaveformTypeSin,
}

def stimx_waveform_factory(conf):
    """
    Factory to create correct waveform class based on conf.
    Supports json config in conf as well as string pointer to a file.
    :rtype: BaseWaveformType
    """
    waveform_conf = conf["extracellular_stimelectrode"]["waveform"]

    if type(waveform_conf) in (str,unicode): # if waveform_conf is str or unicode assume to be name of file in stim_dir
        # waveform_conf = str(waveform_conf)   # make consistent
        wf_dir = conf["manifest"]["$STIM_DIR"] + "/"
        filetype = waveform_conf.split('.')[-1] if '.' in waveform_conf else None

        if filetype == 'csv':
            return WaveformCustom(waveform_conf,wf_dir)
        elif filetype == 'json':
            with open(wf_dir + waveform_conf, 'r') as f:
                waveform_conf = json.load(f)
        else:
            print "Warning: unknwon filetype for waveform"

    shape_key = waveform_conf["shape"].lower()

    if shape_key not in shape_classes:
        print "Warning: waveform shape not known" # throw error?

    Constructor = shape_classes[shape_key]

    return Constructor(waveform_conf)

def iclamp_waveform_factory(conf):
    """
    Factory to create correct waveform class based on conf.
    Supports json config in conf as well as string pointer to a file.
    :rtype: BaseWaveformType
    """
    iclamp_waveform_conf = conf["iclamp"]

    shape_key = iclamp_waveform_conf["shape"].lower()

    if shape_key not in shape_classes:
        print "Warning: iclamp waveform shape not known" # throw error?

    Constructor = shape_classes[shape_key]

    return Constructor(iclamp_waveform_conf)