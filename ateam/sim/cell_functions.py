"""Additional axon processing options for all-active models"""


from neuron import h
from bmtk.simulator.bionet.pyfunction_cache import add_cell_processor

from bmtk.simulator.bionet.io_tools import io
from bmtk.simulator.bionet.default_setters.cell_models import fix_axon_peri, fix_axon_allactive, \
                    set_params_allactive, fix_axon_perisomatic_directed,get_axon_direction
import numpy as np


def fix_axon_allactive_bpopt(hobj):
  """Replace reconstructed axon with a stub which is consistent with BluePyOpt
  Parameters
  ----------
  hobj: instance of a Biophysical template
      NEURON's cell object
  """
  # Set axon diameter to the first section diameter unless reconstructed axon is 60 microns long
  
  nsec = 0
  for sec in hobj.all:
     section_name = sec.name().split(".")[1][:4]
     if section_name == 'axon':
         nsec += 1
  
  if nsec == 0:
      axon_diams = [1,1]
  elif nsec == 1:
      axon_diams = [hobj.axon[0].diam, hobj.axon[0].diam]
  else:
      axon_diams = [hobj.axon[0].diam, hobj.axon[0].diam]
      h.distance(sec=hobj.soma[0])   # need this to set all distances relative to soma (not sure if from center?)
      for sec in hobj.all:
         section_name = sec.name().split(".")[1][:4]
         if section_name == 'axon':
             for seg in sec:
               if h.distance(seg.x) > 60:
                 axon_diams[1] = sec.diam
         

  for sec in hobj.axon:
      h.delete_section(sec=sec)

  h.execute('create axon[2]', hobj)
  for index, sec in enumerate(hobj.axon):
      sec.L = 30
      sec.diam = axon_diams[index]  # 1

      hobj.axonal.append(sec=sec)
      hobj.all.append(sec=sec)  # need to remove this comment

  hobj.axon[0].connect(hobj.soma[0], 1.0, 0)
  hobj.axon[1].connect(hobj.axon[0], 1.0, 0)

  h.define_shape()

def fix_axon_peri_ani(hobj):
    """Replace reconstructed axon with a stub

    :param hobj: hoc object
    """
    for sec in hobj.axon:
        h.delete_section(sec=sec)

    h.execute('create axon[2]', hobj)

    for sec in hobj.axon:
        sec.L = 30
        sec.diam = 1
        hobj.axonal.append(sec=sec)
        hobj.all.append(sec=sec)  # need to remove this comment

    hobj.axon[0].connect(hobj.soma[0], 1, 0)
    hobj.axon[1].connect(hobj.axon[0], 1, 0)

    h.define_shape()

def fix_axon_peri_ani_directed(hobj):
    # io.log_info('Fixing Axon like perisomatic')
    all_sec_names = []
    for sec in hobj.all:
        all_sec_names.append(sec.name().split(".")[1][:4])

    if 'axon' not in all_sec_names:
        io.log_exception('There is no axonal recostruction in swc file.')
    else:
        beg1, end1, beg2, end2 = get_axon_direction(hobj)

    for sec in hobj.axon:
        h.delete_section(sec=sec)
    h.execute('create axon[2]', hobj)

    h.pt3dadd(beg1[0], beg1[1], beg1[2], 1, sec=hobj.axon[0])
    h.pt3dadd(end1[0], end1[1], end1[2], 1, sec=hobj.axon[0])
    hobj.all.append(sec=hobj.axon[0])
    h.pt3dadd(beg2[0], beg2[1], beg2[2], 1, sec=hobj.axon[1])
    h.pt3dadd(end2[0], end2[1], end2[2], 1, sec=hobj.axon[1])
    hobj.all.append(sec=hobj.axon[1])

    hobj.axon[0].connect(hobj.soma[0], 1.0, 0)
    hobj.axon[1].connect(hobj.axon[0], 1.0, 0)

    hobj.axon[0].L = 30.0
    hobj.axon[1].L = 30.0

    h.define_shape()

    for sec in hobj.axon:
        # print "sec.L:", sec.L
        if np.abs(30-sec.L) > 0.0001:
            io.log_exception('Axon stub L is less than 30')



def set_params_allactive_AIS_seg(hobj, params_dict,**kwargs):
    '''
        Sets the properties of section elements defined in kwargs from 
        the Section template ('soma', 'axon', 'apic', 'dend')
        defined in kwargs. If no section template is defined the section elements
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
     
    select_sec_names = kwargs.get('select_section_names', [])
    select_section =  kwargs.get('select_section','')
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
            if sec.name() not in select_sec_names:
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


def add_ais_segment(hobj):
    
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

    hobj.axon[0].connect(hobj.soma[0], 1.0, 0)
    hobj.axon[1].connect(hobj.axon[0], 1.0, 0)
    hobj.axon[2].connect(hobj.axon[1], 1.0, 0)

    h.define_shape()
    return ais_sec_names


def allactive_ais_passive(hobj, cell, dynamics_params):
    # Adds a segment between the AIS and soma,
    # with passive properties from axon
    phantom_sec_names = add_ais_segment(hobj)
    set_params_allactive_AIS_seg(hobj, dynamics_params, \
            select_section_names = phantom_sec_names)
    return hobj

def allactive_ais_somatic(hobj, cell, dynamics_params):
    # Adds a segment between the AIS and soma,
    # with all properties matching the soma
    phantom_sec_names = add_ais_segment(hobj)
    set_params_allactive_AIS_seg(hobj, dynamics_params,\
                select_section_names = phantom_sec_names,
                select_section = 'soma')
    
    return hobj

def aibs_allactive_bpopt_axon(hobj, cell, dynamics_params):
   # Write custom function for replacing axon with stub
   # custom_axons_cut(hobj)
   fix_axon_allactive_bpopt(hobj)
   set_params_allactive(hobj, dynamics_params)
   return hobj


def aibs_allactive_ani(hobj, cell, dynamics_params):
    return aibs_allactive_stub_axon(hobj, cell, dynamics_params)

def aibs_allactive_stub_axon(hobj, cell, dynamics_params):
   fix_axon_peri_ani(hobj) # Replace axon with a stub 60 micron with 1 micron diameter
   set_params_allactive(hobj, dynamics_params)
   return hobj

def aibs_allactive_ani_directed(hobj, cell, dynamics_params):
   fix_axon_peri_ani_directed(hobj) # Replace axon with a stub 60 micron with 1 micron diameter
   set_params_allactive(hobj, dynamics_params)
   return hobj


add_cell_processor(allactive_ais_somatic)
add_cell_processor(allactive_ais_passive)
add_cell_processor(aibs_allactive_bpopt_axon)
add_cell_processor(aibs_allactive_stub_axon)
add_cell_processor(aibs_allactive_ani)
add_cell_processor(aibs_allactive_ani_directed)
