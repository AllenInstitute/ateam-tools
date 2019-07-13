from neuron import h
from bmtk.simulator.bionet.io_tools import io
from bmtk.simulator.bionet.default_setters.cell_models import fix_axon_peri, fix_axon_allactive, \
                    set_params_allactive, fix_axon_perisomatic_directed,get_axon_direction
import numpy as np

def set_params_allactive_AIS_seg(hobj, params_dict, select_section_names=[], select_section=''):
    '''
        Sets the properties of section elements defined in select_section_names from 
        the Section template ('soma', 'axon', 'apic', 'dend')
        defined in select_section. If no section template is defined the section elements
        were only set passive properties
    '''

    passive = params_dict['passive'][0]
    genome = params_dict['genome']
    conditions = params_dict['conditions'][0]

    section_map = {}
    for sec in hobj.all:
        section_name = sec.name().split(".")[1][:4]
        if section_name in section_map:
            section_map[section_name].append(sec)
        else:
            section_map[section_name] = [sec]

    for sec in hobj.all:
        sec.insert('pas')
        # sec.insert('extracellular')

    if 'e_pas' in passive:
        e_pas_val = passive['e_pas']
        for sec in hobj.all:
            for seg in sec:
                seg.pas.e = e_pas_val

    if 'ra' in passive:
        ra_val = passive['ra']
        for sec in hobj.all:
            sec.Ra = ra_val

    if 'cm' in passive:
        for cm_dict in passive['cm']:
            cm = cm_dict['cm']
            for sec in section_map.get(cm_dict['section'], []):
                sec.cm = cm
     
    genome_dict_select_list = []
    sec_select_list = []
    
    for genome_dict in genome:
        g_section = genome_dict['section']
        if genome_dict['section'] == 'glob':
            io.log_warning("There is a section called glob, probably old json file")
            continue

        g_value = float(genome_dict['value'])
        g_name = genome_dict['name']
        g_mechanism = genome_dict.get("mechanism", "")
        for sec in section_map.get(g_section, []):
            if sec.name() not in select_section_names:
                if g_mechanism != "":
                    sec.insert(g_mechanism)
                setattr(sec, g_name, g_value)
            else:
                sec_select_list.append(sec)
                
        
        if select_section != '' and select_section == g_section:
            genome_dict_select_list.append(genome_dict)
        elif select_section == '' and g_mechanism == '' and g_section == 'axon':
            genome_dict_select_list.append(genome_dict)
    
    sec_select_list = list(set(sec_select_list))        
    
    if len(genome_dict_select_list) > 0:    
        for genome_dict_select in genome_dict_select_list:
            g_section = genome_dict_select['section']
            
            g_value = float(genome_dict_select['value'])
            g_name = genome_dict_select['name']
            g_mechanism = genome_dict_select.get("mechanism", "")
            
            for sec in sec_select_list:
                if g_mechanism != "":
                    sec.insert(g_mechanism)
                setattr(sec, g_name, g_value)

        
    for erev in conditions['erev']:
        erev_section = erev['section']
        erev_ena = erev['ena']
        erev_ek = erev['ek']

        if erev_section in section_map:
            for sec in section_map.get(erev_section, []):
                if h.ismembrane('k_ion', sec=sec) == 1:
                    setattr(sec, 'ek', erev_ek)
                if h.ismembrane('na_ion', sec=sec) == 1:
                    setattr(sec, 'ena', erev_ena)
        else:
            io.log_warning("Can't set erev for {}, section array doesn't exist".format(erev_section))

def add_ais_segment(hobj, directed=False):
    if directed:
        try:
            beg, end = get_axon_direction(hobj)[:2]
        except:
            io.log_exception('Error finding axon direction from swc.')

    axon_diams = [1,1,1]
    for sec in hobj.axon:
      h.delete_section(sec=sec)
     
    ais_sec_names = []    
    h.execute('create axon[3]', hobj)
    for index, sec in enumerate(hobj.axon):
        if index == 0:
            #sec.L = 20
            ais_sec_names.append(sec.name())
        #else:
            #sec.L = 30
        sec.L = 30
            
        sec.diam = axon_diams[index]  # 1

        hobj.axonal.append(sec=sec)
        hobj.all.append(sec=sec)  # need to remove this comment
        if directed:
            h.pt3dadd(beg[0], beg[1], beg[2], sec.diam, sec=sec)
            h.pt3dadd(end[0], end[1], end[2], sec.diam, sec=sec)

    hobj.axon[0].connect(hobj.soma[0], 1.0, 0)
    hobj.axon[1].connect(hobj.axon[0], 1.0, 0)
    hobj.axon[2].connect(hobj.axon[1], 1.0, 0)

    h.define_shape()
    return ais_sec_names

def allactive_ais_passive(hobj, cell, dynamics_params):
    # Adds a segment between the AIS and soma,
    # with passive properties from axon
    phantom_sec_names = add_ais_segment(hobj, directed=True)
    set_params_allactive_AIS_seg(hobj, dynamics_params, \
            select_section_names = phantom_sec_names)
    return hobj

def allactive_ais_somatic(hobj, cell, dynamics_params):
    # Adds a segment between the AIS and soma,
    # with all properties matching the soma
    phantom_sec_names = add_ais_segment(hobj, directed=True)
    set_params_allactive_AIS_seg(hobj, dynamics_params,\
                select_section_names = phantom_sec_names,
                select_section = 'soma')
    
    return hobj
