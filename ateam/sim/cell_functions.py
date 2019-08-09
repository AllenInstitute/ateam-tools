"""Additional axon processing options for all-active models"""


from neuron import h
from bmtk.simulator.bionet.pyfunction_cache import add_cell_processor

from bmtk.simulator.bionet.io_tools import io
from bmtk.simulator.bionet.default_setters.cell_models import fix_axon_peri, fix_axon_allactive, \
                    set_params_allactive, fix_axon_perisomatic_directed,get_axon_direction
from .ais_functions import allactive_ais_somatic, allactive_ais_passive
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

def fix_axon_allactive_bpopt_directed(hobj):
    nsec = 0
    for sec in hobj.all:
        section_name = sec.name().split(".")[1][:4]
        if section_name == 'axon':
            nsec += 1
  
    if nsec == 0:
        io.log_exception('There is no axonal recostruction in swc file.')
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
    hobj.axon[0].diam = axon_diams[0]
    hobj.axon[1].diam = axon_diams[1]
    
    # Needed to update section start points after length changed
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

    # Needed to update section start points after length changed
    h.define_shape()

    for sec in hobj.axon:
        # print "sec.L:", sec.L
        if np.abs(30-sec.L) > 0.0001:
            io.log_exception('Axon stub L is less than 30')

def aibs_allactive_bpopt_axon(hobj, cell, dynamics_params):
   # Write custom function for replacing axon with stub
   # custom_axons_cut(hobj)
   fix_axon_allactive_bpopt(hobj)
   set_params_allactive(hobj, dynamics_params)
   return hobj

def aibs_allactive_bpopt_axon_directed(hobj, cell, dynamics_params):
    fix_axon_allactive_bpopt_directed(hobj)
    set_params_allactive(hobj, dynamics_params)
    return hobj



def aibs_allactive_ani(hobj, cell, dynamics_params):
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
add_cell_processor(aibs_allactive_bpopt_axon_directed)
add_cell_processor(aibs_allactive_ani)
add_cell_processor(aibs_allactive_ani_directed)
