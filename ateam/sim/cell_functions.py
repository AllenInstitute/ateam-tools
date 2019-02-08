"""Additional axon processing options for all-active models"""


from neuron import h
from bmtk.simulator.bionet.pyfunction_cache import add_cell_processor
from bmtk.simulator.bionet.default_setters.cell_models import fix_axon_peri, fix_axon_allactive, set_params_allactive

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

def add_ais_segment(hobj, name):
    ais = h.Section(name=name, cell=hobj)
    ais.L = 10
    ais.diam = 1
    hobj.all.append(sec=ais)
    h.disconnect(sec=hobj.axon[0])
    ais.connect(hobj.soma[0], 1.0, 0)
    hobj.axon[0].connect(ais, 1.0, 0)
    h.define_shape()
    return ais

def allactive_ais_passive(hobj, cell, dynamics_params):
    # Adds a segment between the AIS and soma,
    # with passive properties (default cm=1 and Rm and Ra from params)
    fix_axon_peri(hobj)
    ais = add_ais_segment(hobj, name='ais')
    set_params_allactive(hobj, dynamics_params)
    return hobj

def allactive_ais_somatic(hobj, cell, dynamics_params):
    # Adds a segment between the AIS and soma,
    # with active properties matching the soma
    fix_axon_peri(hobj)
    # the soma2 name will be treated as somatic by set_params_allactive
    # (checks only the first 4 chars of names)
    ais = add_ais_segment(hobj, name='soma2')
    set_params_allactive(hobj, dynamics_params)
    return hobj

def aibs_allactive_bpopt_axon(hobj, cell, dynamics_params):
   # Write custom function for replacing axon with stub
   # custom_axons_cut(hobj)
   fix_axon_allactive_bpopt(hobj)
   set_params_allactive(hobj, dynamics_params)
   return hobj

def aibs_allactive(hobj, cell, dynamics_params):
   
   fix_axon_peri(hobj) # Replace axon with a stub 60 micron with 1 micron diameter
   set_params_allactive(hobj, dynamics_params)
   return hobj


add_cell_processor(allactive_ais_somatic)
add_cell_processor(allactive_ais_passive)
add_cell_processor(aibs_allactive_bpopt_axon)
add_cell_processor(aibs_allactive)
