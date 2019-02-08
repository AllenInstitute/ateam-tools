#!/shared/utils.x86_64/python-2.7/bin/python
import os
#import shutil
import sys
import glob
import re
import psycopg2
import psycopg2.extras
import prep_upright
import filecmp
from neuron_morphology import swc

from swank_common import *
import mouse

########################################################################
# database interface code
try:
    conn_string = "host='limsdb2' dbname='lims2' user='atlasreader' password='atlasro'"
    conn = psycopg2.connect(conn_string)
    cursor = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
except:
    print "-------------------------"
    print "Unable to connect to LIMS"
    raise

def specimen_id_from_name(name):
    global cursor
    cursor.execute("SELECT id FROM specimens WHERE name ilike '%%%s%%'" % name)
    result = cursor.fetchall()
    if len(result) > 0 and len(result[0]) > 0:
        return str(result[0][0])

CATEGORY_MOUSE = "Mouse"
CATEGORY_HUMAN = "Human"
CATEGORY_UNKNOWN = "Unrecognized"
CATEGORY_LGN = "LGN"

# returns whether cell is human or mouse, or LGN
def specimen_category(spec_id):
    global cursor
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
    m = re.compile("\d{6}\.\d\d\.\d\d\.\d\d")
    obj = m.search(filename)
    if obj is not None:
        return specimen_id_from_name(obj.group())
    return None

# copy from master pixel storage to upright directory (or non-upright) and
#   convert to micron scale
def copy_to_upright(in_filename, out_filename, spec_id):
    nrn = swc.read_swc(in_filename)
    if len(in_filename.split("_p")) > 1:
        # convert to micron scale
        for n in nrn.node_list:
            n.x *= 0.113
            n.y *= 0.113
            n.z *= 3 * 0.113
            n.radius *= 0.113
    elif len(in_filename.split("_m")) > 1:
        # scale z axis
        for n in nrn.node_list:
            n.z *= 3
    else:
        print "file %s not a recognized pixel or micron file" % in_filename
    # set soma to origin
    soma = nrn.soma_root()
    dx = soma.x
    dy = soma.y
    dz = soma.z
    for n in nrn.node_list:
        n.x -= dx
        n.y -= dy
        n.z -= dz
    # apply affine
    try:
        aff = prep_upright.calculate_transform(cursor, spec_id, nrn)
    except:
        aff = None
    # see where cell should be routed to
    target = specimen_category(spec_id)
    if target == CATEGORY_MOUSE:
        if aff is not None:
            print("Saving %s to upright directory" % out_filename)
            nrn.apply_affine(aff)
            master_file = UPRIGHT_DIR + out_filename
            nrn.save(master_file)
            mouse_file = MOUSE_INBOUND_DIR + out_filename
            nrn.save(mouse_file)
        else:
            print("Saving %s to non-upright directory" % out_filename)
            nrn.save(RAW_DIR + out_filename)
            nrn.save(MOUSE_INBOUND_NON_UPRIGHT_DIR + out_filename)
    elif target == CATEGORY_HUMAN:
        print "Human cell - doing nothing at present"
    elif target == CATEGORY_LGN:
        print "LGN cell - doing nothing at present"
    else:
        print "%s is from unrecognized specimen category" % spec_id

# perform upright on cells and copy to appropriate target project,
#   then delete from upload directory
swcs = glob.glob(UPLOAD_DIR + "*.swc")
print("Found %d files in upload directory" % len(swcs))
for s in swcs:
    spec_id = get_spec_id(s)
    print "Preprocessing %s:    %s" % (spec_id, s)
    try:
        # copy pixel file to
        filename = generate_outfile_name(spec_id)
        # create uprighted micron version
        copy_to_upright(s, filename, spec_id)
        # when process complete, delete input file
        os.remove(s)
    except:
        print("Error copying/deleting '%s'" % s)
        raise

print("")

# process mouse data
mouse.process_mouse_data()
# TODO process human data
# TODO process lgn data

