import neurom as nm
import numpy as np
from neurom import morphmath as mm
from neurom.fst.sectionfunc import section_path_length
from neurom.core import Tree, iter_neurites, iter_sections, iter_segments, NeuriteType
from neurom.core.types import tree_type_checker

# Neurite types of interest
NEURITES = (nm.NeuriteType.all,
             nm.NeuriteType.apical_dendrite,
             nm.NeuriteType.basal_dendrite,)

def process_morph_file(morph_file):
    nrn = nm.load_neuron(morph_file)

    morph_stats = {}

    morph_stats['soma_suface'] = nm.get('soma_surface_areas', nrn)[0]
    morph_stats['soma_radius'] = np.mean(nm.get('soma_radii', nrn))

    neurite_filter = tree_type_checker(nm.NeuriteType.apical_dendrite, nm.NeuriteType.basal_dendrite)
    morph_stats['length'] = sum(mm.segment_length(s) for s in nm.iter_segments(nrn, neurite_filter=neurite_filter))
    morph_stats['area'] = sum(mm.segment_area(s) for s in nm.iter_segments(nrn, neurite_filter=neurite_filter))
    morph_stats['volume'] = sum(mm.segment_volume(s) for s in nm.iter_segments(nrn, neurite_filter=neurite_filter))
    morph_stats['max_eff_length'] = max_fcn_terminals(section_eff_path_length, nrn, neurite_filter=neurite_filter)
    morph_stats['max_path_length'] = max_fcn_terminals(section_path_length, nrn, neurite_filter=neurite_filter)
    # Morph stats
    # for nrn_type_ in NEURITES:
    #     type_name = str(nrn_type_).split('.')[1]
    #     morph_stats['length_'+type_name] = np.sum(nm.get('segment_lengths', nrn, neurite_type=nrn_type_))
    #     morph_stats['area_'+type_name] = sum(mm.segment_area(s) for s in nm.iter_segments(nrn, neurite_filter=tree_type_checker(nrn_type_)))
    #     morph_stats['volume_'+type_name] = sum(mm.segment_volume(s) for s in nm.iter_segments(nrn, neurite_filter=tree_type_checker(nrn_type_)))
    #     morph_stats['taper_rate_'+type_name] = np.mean([mm.segment_taper_rate(s) for s in nm.iter_segments(
    #             nrn, neurite_filter=tree_type_checker(nrn_type_))])
    #     morph_stats['max_eff_length_'+type_name] = max_fcn_terminals(section_eff_path_length, nrn, neurite_type=nrn_type_)
    #     morph_stats['max_path_length_'+type_name] = max_fcn_terminals(section_path_length, nrn, neurite_type=nrn_type_)
    return morph_stats

def max_fcn_terminals(sec_fcn, neurites, neurite_filter):
    '''Get the path lengths to each terminal point per neurite in a collection'''
    return max(sec_fcn(s)
            for n in iter_neurites(neurites, filt=neurite_filter)
            for s in iter_sections(n, iterator_type=Tree.ileaf))

def section_eff_path_length(section):
    '''Path length from section to root'''
    return sum(s.length / np.sqrt(section_radius(s)) for s in section.iupstream())

def section_radius(s):
    radius = np.mean([mm.segment_radius(seg) for seg in iter_segments(s)])
    return radius
