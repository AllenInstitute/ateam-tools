import logging
import os
import psycopg2
import psycopg2.extras
import prep_upright
from contextlib import contextmanager
import re
import numpy as np

CATEGORY_MOUSE = "Mouse"
CATEGORY_HUMAN = "Human"
CATEGORY_UNKNOWN = "Unrecognized"
CATEGORY_LGN = "LGN"

@contextmanager
def open_lims():
    conn_string = "host='limsdb2' dbname='lims2' user='atlasreader' password='atlasro'"
    conn = psycopg2.connect(conn_string)
    cursor = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

    yield cursor

    cursor.close()

def specimen_id_from_image_series_id(iser_id):
    with open_lims() as cursor:
        cursor.execute("SELECT specimen_id from image_series where id = %s" % str(iser_id))
        result = cursor.fetchall()
        if len(result) > 0 and len(result[0]) > 0:
            return str(result[0][0])

def specimen_id_from_name(name):
    with open_lims() as cursor:
        cursor.execute("SELECT id FROM specimens WHERE name ilike '%%%s%%'" % name)
        result = cursor.fetchall()
        if len(result) > 0 and len(result[0]) > 0:
            return str(result[0][0])

def generate_outfile_name(spec_id):
    return str(spec_id) + ".swc"

# returns whether cell is human or mouse, or LGN
def specimen_category(spec_id):
    sql = """
    select o.name from specimens cell
    join donors d on d.id = cell.donor_id
    join organisms o on o.id = d.organism_id
    where cell.id =
    """
    lgn_sql = """
    select o.name from specimens cell
    join specimens slice on slice.id = cell.parent_id
    join donors d on d.id = cell.donor_id
    join organisms o on o.id = d.organism_id
    join tissue_processings tp on tp.id = slice.tissue_processing_id
    join tcp_sops tsop on tsop.id = tp.tcp_sop_id
    where tsop.name ilike '%LGN%'
    and cell.id =
    """

    with open_lims() as cursor:
        cursor.execute(lgn_sql + spec_id)
        result = cursor.fetchall()
        if len(result) > 0 and len(result[0]) > 0:
            return CATEGORY_LGN
        #
        cursor.execute(sql + spec_id)
        result = cursor.fetchall()
        
    if len(result) > 0 and len(result[0]) > 0:
        species = result[0][0]
        if species == "Mus musculus":
            return CATEGORY_MOUSE
        elif species == "Homo Sapiens":
            return CATEGORY_HUMAN
    return CATEGORY_UNKNOWN


#
########################################################################

def get_spec_id(filename):
    m = re.compile("\d{9}\.swc")
    obj = m.search(filename)
    if obj is not None:
        try:
            return str(obj.group()[:-4])
        except:
            pass

    m = re.compile("\d{9}")
    obj = m.search(filename)
    if obj is not None:
        return specimen_id_from_image_series_id(obj.group())
    return None


def align_and_scale_neuron(spec_id, nrn, units):
    if units == 'pixels':
        PIXEL_SCALE = dict(x=0.1144, 
                           y=0.1144, 
                           z=0.28, 
                           r=0.1144)
        logging.debug("scaling by %s" % str(PIXEL_SCALE))

        # convert to micron scale
        for n in nrn.node_list:
            n.x *= PIXEL_SCALE['x']
            n.y *= PIXEL_SCALE['y']
            n.z *= PIXEL_SCALE['z']
            n.radius *= PIXEL_SCALE['r']

    elif units == 'microns':
        MICRON_Z_SCALE = 1.0
        logging.debug("scaling z axis by %f" % MICRON_Z_SCALE)
        # scale z axis
        for n in nrn.node_list:
            n.z *= MICRON_Z_SCALE
    else:
        logging.warning("neuron has unknown units (%s) - no scaling applied." % units)

    # set soma to origin
    soma = nrn.soma_root()
    if soma is None:
        logging.warning("could not find soma, using root")
        soma = nrn.root()
    dx = soma.x
    dy = soma.y
    dz = soma.z
    for n in nrn.node_list:
        n.x -= dx
        n.y -= dy
        n.z -= dz
    # apply affine
    try:
        with open_lims() as cursor:
            aff = prep_upright.calculate_transform(cursor, spec_id, nrn)

    except Exception as e:
        logging.error("no upright transform available for specimen %s" % spec_id)
        logging.warning("no upright transform available for specimen %s" % spec_id)
        aff = None
    
    
    logging.warning("upright transform result %s" % spec_id)
    return aff

# copy from master pixel storage to upright directory (or non-upright) and
#   convert to micron scale
    

def copy_to_upright(nrn, aff, base_file_name, 
                    upright_dir, upright_to_process_dir,
                    non_upright_dir, non_upright_to_process_dir):
    # see where cell should be routed to
    if aff is not None:
        nrn.apply_affine(aff)
        upright_file = os.path.join(upright_dir, base_file_name)
        logging.debug("Saving %s to upright directory" % upright_file)
        nrn.save(upright_file)
        inbound_file = os.path.join(upright_to_process_dir, base_file_name)
        logging.debug("Saving %s to upright to-process directory" % inbound_file)
        nrn.save(inbound_file)
    else:
        raw_file = os.path.join(non_upright_dir, base_file_name)
        logging.debug("Saving %s to non-upright directory" % raw_file)
        nrn.save(raw_file)
        inbound_file = os.path.join(non_upright_to_process_dir, base_file_name)
        logging.debug("Saving %s to non-upright to-process directory" % raw_file)
        nrn.save(inbound_file)


def swc_units(file_name):
    if len(file_name.split("_p")) > 1:
        return "pixels"
    elif len(file_name.split("_m")) > 1:
        return "microns"
    else:
        logging.error("file %s not a recognized pixel or micron file" % file_name)

def convert_from_titan_linux(file_name):
    # Lookup table mapping project to program
    project_to_program= {
        "neuralcoding": "braintv", 
        '0378': "celltypes",
        'conn': "celltypes",
        'ctyconn': "celltypes",
        'humancelltypes': "celltypes",
        'mousecelltypes': "celltypes",
        'shotconn': "celltypes",
        'synapticphys': "celltypes",
        'whbi': "celltypes",
        'wijem': "celltypes"
    }
    # Tough intermediary state where we have old paths
    # being translated to new paths
    m = re.match('/projects/([^/]+)/vol1/(.*)', file_name)
    if m:
        newpath = os.path.normpath(os.path.join(
            '/allen',
            'programs',
            project_to_program.get(m.group(1),'undefined'),
            'production',
            m.group(1),
            m.group(2)
        ))
        return newpath
    return file_name

def storage_paths(storage_dir):
    morph_base = os.path.join(storage_dir, 'morphology_data')

    return {
        CATEGORY_MOUSE: {
            'autotrace': os.path.join(morph_base, 'mouse', 'autotrace'),
            'upright': os.path.join(morph_base, 'mouse', 'swc'),
            'non_upright': os.path.join(morph_base, 'mouse', 'swc', 'non_upright'),
            'upright_to_process': os.path.join(morph_base, 'mouse', 'to_process'),
            'non_upright_to_process': os.path.join(morph_base, 'mouse', 'to_process', 'non_upright'), 
            'features': os.path.join(morph_base, 'mouse', 'features'),
            'html': os.path.join(morph_base, 'mouse', 'html')
            },
        CATEGORY_HUMAN: {
            'autotrace': os.path.join(morph_base, 'human', 'autotrace'),
            'upright': os.path.join(morph_base, 'human', 'swc'),
            'non_upright': os.path.join(morph_base, 'human', 'swc', 'non_upright'),
            'upright_to_process': os.path.join(morph_base, 'human', 'to_process'),
            'non_upright_to_process': os.path.join(morph_base, 'human', 'to_process', 'non_upright'),  
            'features': os.path.join(morph_base, 'human', 'features'),
            'html': os.path.join(morph_base, 'human', 'html')
            },
        CATEGORY_LGN: {
            'autotrace': os.path.join(morph_base, 'lgn', 'autotrace'),
            'upright': os.path.join(morph_base, 'lgn', 'swc'),
            'non_upright': os.path.join(morph_base, 'lgn', 'swc', 'non_upright'),
            'upright_to_process': os.path.join(morph_base, 'lgn', 'to_process'),
            'non_upright_to_process': os.path.join(morph_base, 'lgn', 'to_process', 'non_upright'),  
            'features': os.path.join(morph_base, 'lgn', 'features'),
            'html': os.path.join(morph_base, 'lgn', 'html')
            },
        }

def autotrace_path(spec_id):
    with open_lims() as cursor:
        query = """
            SELECT  (cell.storage_directory || wkf.filename) 
            FROM neuron_reconstructions nr
            JOIN well_known_files wkf 
                ON wkf.attachable_id = nr.id
            JOIN specimens cell 
                ON cell.id = nr.specimen_id
            WHERE wkf.filename ILIKE '%_m.swc'
                AND nr.manual = false
                AND wkf.filename NOT ILIKE '%marker%'
                AND cell.id=
        """
        cursor.execute(query + spec_id)
        result = cursor.fetchall()

    if len(result) > 0 and len(result[0]) > 0:
        return convert_from_titan_linux(str(result[0][0]))
    else:
        raise IOError("Unable to fetch autotrace for " + spec_id)
