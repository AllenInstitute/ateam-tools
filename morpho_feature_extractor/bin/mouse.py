#!/shared/utils.x86_64/python-2.7/bin/python
import logging
import os
import shutil
import sys
import glob
import psycopg2
import psycopg2.extras
import time
import filecmp
from neuron_morphology import swc
from neuron_morphology.features.feature_extractor import MorphologyFeatures
import neuron_morphology.visualization.morphvis as morphvis
import swank
import multiprocessing as mp


def specimen_id_from_name(name):
    with swank.open_lims() as cursor:
        cursor.execute("SELECT id FROM specimens WHERE name ilike '%%%s%%'" % name)
        result = cursor.fetchall()

    if len(result) > 0 and len(result[0]) > 0:
        return str(result[0][0])

    raise Exception("Can't find specimen ID for name including '%s'" % name)

def fetch_specimen_record(spec_id):
    with swank.open_lims() as cursor:
        base_sql = ""
        base_sql += "with dendrite_type as  \n"
        base_sql += "( \n"
        base_sql += "  select sts.specimen_id, st.name  \n"
        base_sql += "  from specimen_tags_specimens sts \n"
        base_sql += "  join specimen_tags st on sts.specimen_tag_id = st.id \n"
        base_sql += "  where st.name like 'dendrite type%s' \n"
        base_sql += ") \n"
        #
        name_sql = base_sql
        name_sql += "SELECT spec.id, spec.name, dt.name, str.name, wkf.filename, wkf.storage_directory \n"
        name_sql += "FROM specimens spec \n"
        name_sql += "LEFT JOIN structures str on spec.structure_id = str.id \n"
        name_sql += "LEFT JOIN dendrite_type dt on dt.specimen_id = spec.id \n"
        name_sql += "JOIN neuron_reconstructions nr ON nr.specimen_id=spec.id \n"
        name_sql += "  AND nr.superseded = 'f' AND nr.manual = 't' \n"
        name_sql += "JOIN well_known_files wkf ON wkf.attachable_id=nr.id \n"
        name_sql += "  AND wkf.attachable_type = 'NeuronReconstruction' \n"
        name_sql += "JOIN cell_soma_locations csl ON csl.specimen_id=spec.id \n"
        name_sql += "JOIN well_known_file_types wkft \n"
        name_sql += "  ON wkft.id=wkf.well_known_file_type_id \n"
        name_sql += "WHERE spec.id="
        name2_sql = " AND wkft.name = '3DNeuronReconstruction'; \n"
        #
        cursor.execute(name_sql + spec_id + name2_sql)
        result = cursor.fetchall()

    record = None
    if len(result) > 0:
        record = {}
        record["spec_id"] = result[0][0]
        record["spec_name"] = result[0][1]
        if result[0][2] is not None:
            record["dend_type"] = result[0][2].replace(","," ")
        else:
            record["dend_type"] = ""
        if result[0][3] is not None:
            record["location"] = result[0][3].replace(","," ")
        else:
            record["location"] = ""
        record["filename"] = result[0][4]
        record["path"] = result[0][5]
        record["cre_state"] = fetch_reporter_state(spec_id)
    return record

def fetch_reporter_state(spec_id):
    with swank.open_lims() as cursor:
        cre_sql = """
            WITH drivers(donor_id, name) as
            (
                with drivers_table(donor_id, name) as (
                select distinct d.id, g.name from donors d
                join donors_genotypes d2g on d2g.donor_id = d.id
                join genotypes g on g.id = d2g.genotype_id
                where g.genotype_type_id = 177835595 --
                order by d.id, g.name )
                select distinct d.id, string_agg(dt.name, ' _AND_ ') from donors d
                join drivers_table dt on dt.donor_id = d.id
                group by d.id
            ),
            reporters(donor_id, name) as
            (
                with reporters_table(donor_id, name) as (
                select distinct d.id, g.name from donors d
                join donors_genotypes d2g on d2g.donor_id = d.id
                join genotypes g on g.id = d2g.genotype_id
                where g.genotype_type_id = 177835597 --
                order by d.id, g.name )
                select distinct d.id, string_agg(rt.name, ' _AND_ ') from donors d
                join reporters_table rt on rt.donor_id = d.id
                group by d.id
            )
            select cr.name from specimens cell
            left join donors d on d.id = cell.donor_id
            left join drivers dr on dr.donor_id = d.id
            left join reporters rep on rep.donor_id = d.id
            left join cell_reporters cr on cr.id = cell.cell_reporter_id
            where cell.id = %d
        """
        cursor.execute(cre_sql % int(spec_id))
        result = cursor.fetchall()

    record = []
    if len(result) == 0:
        reporter = '[NA]'
    else:
        if result[0][0] == 'cre reporter positive':
            reporter = "+"
        elif result[0][0] == 'cre reporter negative':
            reporter = "-"
        elif result[0][0] == 'cre reporter indeterminate':
            reporter = "?"
        elif result[0][0] == 'not available':
            reporter = "(na)"
        else:
            reporter = "(err)"
    return reporter

def fetch_normalized_depth(spec_id):
    with swank.open_lims() as cursor:
        depth_sql = "SELECT normalized_depth FROM cell_soma_locations WHERE specimen_id=%s"
        cursor.execute(depth_sql % spec_id)
        result = cursor.fetchall()

    if len(result) > 0 and len(result[0]) > 0:
        return result[0][0]
    logging.warning("Cannot get normalized depth for %s. Using 0" % spec_id)
    return 0

def cell_name_from_id(spec_id):
    with swank.open_lims() as cursor:
        depth_sql = "SELECT name FROM specimens WHERE id=%s"
        cursor.execute(depth_sql % spec_id)
        result = cursor.fetchall()
    if len(result) > 0 and len(result[0]) > 0:
        return result[0][0]
    return ""


# database interface code
########################################################################
########################################################################
# feature code

def compute_features(swc_file, spec_id, features_dir, swc_dir):
    morph = swc.read_swc(swc_file)
    depth = fetch_normalized_depth(spec_id)
    features = MorphologyFeatures(morph, depth)
# reoganize data
    data = {}
    types = ["axon", "cloud", "dendrite", "basal_dendrite", "apical_dendrite", "all_neurites"]
    data["axon"] = features.axon
    data["cloud"] = features.axon_cloud
    data["dendrite"] = features.dendrite
    data["basal_dendrite"] = features.basal_dendrite
    data["apical_dendrite"] = features.apical_dendrite
    data["all_neurites"] = features.all_neurites
# write csv
# convert .swc to .csv
    outfile = os.path.normpath(os.path.join(features_dir, spec_id + ".txt"))
    logging.debug("Saving features to '%s'" % outfile)
    with open(outfile, 'w') as f:
        record = fetch_specimen_record(spec_id)
        if record is not None:
            logging.debug("found record for %s: %s" % (spec_id, str(record)))
            f.write("%s\t%s\n" % ("spec_name", record["spec_name"]))
            f.write("%s\t%s\n" % ("spec_id", record["spec_id"]))
            f.write("%s\t%s\n" % ("cre_state", record["cre_state"]))
            f.write("%s\t%s\n" % ("dend_type", record["dend_type"]))
            f.write("%s\t%s\n" % ("location", record["location"]))
        else:
            logging.error("no specimen record found for %s" % (spec_id))
        f.write("%s\t%s\n" % ("path", os.path.normpath(os.path.join(swc_dir, spec_id + ".swc"))))

        for t in types:
            f.write("%s\n" % t)
            cols = []
            dataset = data[t]
            for k in dataset:
                cols.append(k)
            cols.sort()
            for c in cols:
                f.write("%s\t%s\n" % (c, dataset[c]))

def consolidate_features_from_ids(spec_ids, csv_name, features_dir):
    with open(os.path.normpath(os.path.join(features_dir, csv_name)), "w") as output:
        header = None
        for spec_id in spec_ids:
            feat = os.path.normpath(os.path.join(features_dir, spec_id + ".txt"))
            # write header line, if first file
            if header is None:
                prefix = ""
                with open(feat, "r") as input:
                    header = ""
                    lines = input.readlines()
                    for line in lines:
                        line = line.strip()
                        toks = line.split('\t')
                        if len(toks) == 1:
                            prefix = toks[0] + "_"
                        else:
                            header += "%s%s," % (prefix, toks[0])
                output.write(header + "\n")
            # write content
            with open(feat, "r") as input:
                lines = input.readlines()
                for line in lines:
                    line = line.strip()
                    toks = line.split('\t')
                    if len(toks) > 1:
                        output.write("%s," % toks[1])
            output.write("\n")

# combine content of all feature files into one csv
def consolidate_features(upright_dir, non_upright_dir, features_dir):
    upright_list = []
    nonupright_list = []

    # make list of uprighted SWCs
    swcs = glob.glob(os.path.normpath(os.path.join(upright_dir, "*.swc")))

    for s in swcs:
        spec_id = swank.get_spec_id(s)
        upright_list.append(spec_id)

    consolidate_features_from_ids(upright_list, "upright_features.csv", features_dir)

    swcs = glob.glob(os.path.normpath(os.path.join(non_upright_dir, "*.swc")))

    for s in swcs:
        spec_id = swank.get_spec_id(s)
        nonupright_list.append(spec_id)

    consolidate_features_from_ids(nonupright_list, "non_upright_features.csv", features_dir)


def run_all_features(upright_dir, non_upright_dir, features_dir):
    swcs = glob.glob(os.path.normpath(os.path.join(upright_dir, "*.swc")))
    logging.debug("Extracting features from %d upright file(s)" % len(swcs))
    for s in swcs:
        spec_id = swank.get_spec_id(s)
        # run features on file
        compute_features(s, spec_id)
    if len(swcs) > 0:
        consolidate_features(upright_dir, non_upright_dir, features_dir)

# feature code
########################################################################
########################################################################
# thumbnail code

PROFILE_THUMB = "_profile.png"
AUTO_THUMB = "_auto.png"
FULL_THUMB = "_full.png"

# creates a thumbnail for the specified SWC file, both from the front 
#   (XY plane) and from the side (ZY plane)
# if the cells is 'raw' (eg, a pixel file), the is_raw flag should be set
def create_profile_view(swc_file, thumb_name, is_raw=False, xy_params=None, zy_params=None, strict_validation=False):
    if xy_params is None:
        xy_params = {}
    if zy_params is None:
        zy_params = {}

    col = morphvis.MorphologyColors()
    nrn = swc.read_swc(swc_file, strict_validation=strict_validation)
    if is_raw:
        soma = nrn.soma_root()
        if soma is None:
            logging.warning("Cannot find soma. Usng root instead.")
            soma = nrn.root()
        dx = soma.x
        dy = soma.y
        dz = soma.z
        for n in nrn.node_list:
            n.x -= dx
            n.y -= dy
            n.z -= dz
            n.z *= 3.0

    # generate thumb with profile view
    WIDTH = 100
    HEIGHT = 120

    img = morphvis.create_image(2*WIDTH, HEIGHT, color=(255, 255, 255))

    thumb = morphvis.create_image(WIDTH, HEIGHT, alpha=True, color=(255, 255, 255, 255))
    xy_params = morphvis.draw_morphology(thumb, nrn, scale_to_fit=not bool(xy_params), colors=col, **xy_params)
    img.paste(thumb, (0, 0), thumb)

    nrn.rotate(90.0)
    thumb = morphvis.create_image(WIDTH, HEIGHT, alpha=True, color=(255, 255, 255, 255))
    zy_params = morphvis.draw_morphology(thumb, nrn, scale_to_fit=not bool(zy_params), colors=col, **zy_params)
    img.paste(thumb, (WIDTH, 0), thumb)
    img.save(thumb_name)

    return xy_params, zy_params

# create all thumbnails for specimen
def generate_thumbnails(spec_id, src_dir, thumb_dir, auto_dir):
    col = morphvis.MorphologyColors()

    # generate thumbs with profile view
    swc_file = os.path.normpath(os.path.join(src_dir, spec_id + ".swc"))
    xy_params, zy_params = create_profile_view(swc_file, os.path.normpath(os.path.join(thumb_dir, spec_id + PROFILE_THUMB)))

    # autotrace profile
    swc_file = os.path.normpath(os.path.join(auto_dir, spec_id + ".swc"))
    if os.path.exists(swc_file):
        create_profile_view(swc_file, os.path.normpath(os.path.join(thumb_dir, spec_id + AUTO_THUMB)),
                            xy_params=xy_params, zy_params=zy_params,
                            strict_validation=None)

    # generate unscaled thumb
    swc_file = os.path.normpath(os.path.join(src_dir, spec_id + ".swc"))
    nrn = swc.read_swc(swc_file)
    dims, low, high = nrn.get_dimensions()
    HEIGHT = 1024
    WIDTH = 1024
    SCALE = 0.8
    img = morphvis.create_image(WIDTH, HEIGHT, color=(255, 255, 255))
    morphvis.draw_morphology(img, nrn, 
                             inset_left=-SCALE*low[0], 
                             inset_top=SCALE*high[1], 
                             scale_factor=SCALE, colors=col)
    #morphvis.draw_morphology(img, nrn, scale_to_fit=True, colors=col)
    img.save(os.path.normpath(os.path.join(thumb_dir, spec_id + FULL_THUMB)))

# thumbnail code
########################################################################
########################################################################
# HTML code

def generate_html(spec_id, upright, html_dir, features_dir):
    with open(os.path.normpath(os.path.join(html_dir, spec_id + ".html")), "w") as html:
        html.write("<html>\n")
        html.write("  <head>\n")
        html.write("    <title>Mouse morphology data</title>\n")
        html.write("  </head>\n")
        html.write("  <body>\n")
        # show thumbnails
        prof_thumb = spec_id + PROFILE_THUMB
        auto_thumb = spec_id + AUTO_THUMB
        full_thumb = spec_id + FULL_THUMB
        html.write('<img src="%s" alt="Profile for %s"><br>\n' % (prof_thumb , spec_id))
        html.write('Profile view<br><hr/>\n')
        html.write('<img src="%s" alt="Autotrace comparison for %s"><br>\n' % (auto_thumb , spec_id))
        html.write('Autotrace profile<br><hr/>\n')
        html.write('<img src="%s" alt="Thumbnail for %s"><br>\n' % (full_thumb , spec_id))
        html.write('Scaled thumbnail<br><hr/>\n')
        # write features
        if upright:
            html.write("<h3>Uprighted cell</h3>\n")
        else:
            html.write("<h3>Non-uprighted cell</h3>\n")
        featfile = os.path.normpath(os.path.join(features_dir, spec_id + ".txt"))
        with open(featfile, 'r') as feat:
            content = feat.readlines()
            for line in content:
                kv = line.split('\t')
                if len(kv) == 1:
                    html.write("<h2>%s</h2>\n" % kv[0].strip())
                else:
                    html.write("&nbsp;%s: %s<br>\n" % (kv[0], kv[1].strip()))
        html.write("  </body>\n")
        html.write("</html\n>")


def generate_index_html(html_dir, upright):
    cell_files = glob.glob(os.path.normpath(os.path.join(html_dir, "[0-9]*.html")))
    with open(os.path.normpath(os.path.join(html_dir, "index.html")), "w") as f:
        f.write("<html>\n")
        f.write("  <head>\n")
        f.write("    <title>Mouse morphology data</title>\n")
        f.write("  </head>\n")
        f.write("  <body>\n")
        # write cell names, ordered by name, grouped by upright-ness
        name_id_map = {}
        name_list = []
        for fname in cell_files:
            spec_id = fname[-14:-5]
            cell_name = cell_name_from_id(spec_id)
            name_id_map[cell_name] = spec_id
            name_list.append(cell_name)
        name_list.sort()
        f.write("<h2>Uprighted cells</h2>\n")
        for name in name_list:
            spec_id = name_id_map[name]
            if spec_id in upright:
                f.write('    <a href="%s.html">%s</a><br>\n' % (spec_id, name))
        f.write("<h2>Non-uprighted cells</h2>\n")
        for name in name_list:
            spec_id = name_id_map[name]
            if spec_id not in upright:
                f.write('    <a href="%s.html">%s</a><br>\n' % (spec_id, name))
        f.write("  </body>\n")
        f.write("</html\n>")

# HTML code
########################################################################
def process_mouse_cell(s, upright,
                       upright_to_process_dir, upright_dir,
                       non_upright_dir, non_upright_to_process_dir,
                       auto_dir, features_dir, html_dir):
    try:
        spec_id = swank.get_spec_id(s)
    except Exception as e:
        print e
        logging.error("Error figuring out specimen ID for file '%s': %s" % (s, e.message))
        return

    logging.info("Processing %s: %s" % (spec_id, s))
    if spec_id in upright:
        in_dir = upright_to_process_dir
        swc_dir = upright_dir
    else:
        in_dir = non_upright_to_process_dir
        swc_dir = non_upright_dir

    try:
        # regenerate filename, w/o path
        filename = swank.generate_outfile_name(spec_id)
        # run features on file
        compute_features(s, spec_id, features_dir, swc_dir)
        # generate figures & html

        generate_thumbnails(spec_id, swc_dir, html_dir, auto_dir)
        generate_html(spec_id, spec_id in upright, html_dir, features_dir)

        # when process complete, delete input file
        os.remove(s)
    except Exception as e:            
        logging.error("Error processing '%s': %s", s, e.message)
        import traceback
        traceback.print_exc()
        raise
        #logging.error("Continuing...")
    
def process_mouse_data(upright_to_process_dir, upright_dir,
                       non_upright_dir, non_upright_to_process_dir,
                       auto_dir, features_dir, html_dir,
                       nprocs=4):
    
    pool = mp.Pool(nprocs)

    upright = {}

    # make master list of upright files
    up = glob.glob(os.path.normpath(os.path.join(upright_to_process_dir, "*.swc")))

    for s in up:
        spec_id = swank.get_spec_id(s)
        upright[spec_id] = s

    swcs = glob.glob(os.path.normpath(os.path.join(upright_dir, "*.swc")))
    for s in swcs:
        spec_id = swank.get_spec_id(s)
        upright[spec_id] = s

    # make list of updated mouse files 
    # files should be uprighted (or not), in microns, and named NNNNNNNNN.swc
    # keep dictionary of upright and nonuprighted cells
    swcs = []
    for s in up:
        spec_id = swank.get_spec_id(s)
        swcs.append(s)
    nonup = glob.glob(os.path.normpath(os.path.join(non_upright_to_process_dir, "*.swc")))
    for s in nonup:
        spec_id = swank.get_spec_id(s)
        swcs.append(s)
    swcs.sort()

    logging.info("Processing %d file(s)" % len(swcs))

    results = []
    for s in swcs:
        result = pool.apply_async(process_mouse_cell, 
                                  args=(s, upright, upright_to_process_dir, upright_dir,
                                        non_upright_dir, non_upright_to_process_dir,
                                        auto_dir, features_dir, html_dir))
        results.append(result)

    pool.close()

    for result in results:
        try:
            result.get()
        except Exception as e:
            logging.exception("Error processing file: %s", e.message)

    pool.join()

    if len(swcs) > 0:
        # consolidate features
        # NOTE: if feature algorithm changes then features need to be
        #   completely rerun for all cells
        consolidate_features(upright_dir, non_upright_dir, features_dir)

        # regenerate index.html
        generate_index_html(html_dir, upright)

if __name__ == "__main__":
    process_mouse_data()

