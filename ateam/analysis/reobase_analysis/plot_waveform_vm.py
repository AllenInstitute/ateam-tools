######################################################
# Authors: Fahimeh Baftizadeh, Taylor Connington
# Date created: 9/1/2017
######################################################

import reobase_analysis.tchelpers as tc
import reobase_analysis.reobase_utils as ru
from reobase_analysis.stimxwaveform import stimx_waveform_factory


# el = 409
# amp = 0.050
# gid = 314900022

def plot(el, amp, cell_gid, model_type):

    example_dir = ru.get_dc_output_dir(cell_gid,el,amp, model_type)
    conf = ru.get_json_from_file(ru.get_config_resolved_path(example_dir, el, amp))
    waveform = stimx_waveform_factory(conf)
    
    title = 'Example $V_m$ and $I_{stim}$ for' + ' el {}, cell {}'.format(el, cell_gid)
    tc.plot_waveform_vm(example_dir, waveform, cell=0, 
                     title=title)
    

# plot(el,amp,gid)