#!/shared/utils.x86_64/python-2.7/bin/python
import logging
import os
#import shutil
import sys
import glob
from neuron_morphology import swc
from neuron_morphology.validation.result import InvalidMorphology
import argparse
import multiprocessing as mp
import swank
import mouse

def upload_swc(s, paths):
    spec_id = swank.get_spec_id(s)

    if spec_id is None:
        raise Exception("could not find specimen id from file name: %s" % s)

    logging.info("Preprocessing %s: %s" % (spec_id, s))


    out_file_name = swank.generate_outfile_name(spec_id)
    units = swank.swc_units(s)
    nrn = swc.read_swc(s)
    aff = swank.align_and_scale_neuron(spec_id, nrn, units)
    category = swank.specimen_category(spec_id)


    autotrace_nrn = None
    try:
        logging.debug("Downloading autotrace path")
        autotrace_file_name = swank.autotrace_path(spec_id)
        autotrace_nrn = swc.read_swc(autotrace_file_name, strict_validation=None)
    except IOError as e:
        logging.error("Autotrace for for specimen %s does not exist: %s", spec_id, e.message)


    if category in paths:
        swank.copy_to_upright(nrn, aff, "%s.swc" % spec_id,
                              upright_dir=paths[category]['upright'],
                              upright_to_process_dir=paths[category]['upright_to_process'],
                              non_upright_dir=paths[category]['non_upright'],
                              non_upright_to_process_dir=paths[category]['non_upright_to_process'])

        if autotrace_nrn is not None:
            swank.align_and_scale_neuron(spec_id, autotrace_nrn, 'microns')
            swank.copy_to_upright(autotrace_nrn, aff, "%s.swc" % spec_id,
                                  upright_dir=paths[category]['autotrace'],
                                  upright_to_process_dir=paths[category]['autotrace'],
                                  non_upright_dir=paths[category]['autotrace'],
                                  non_upright_to_process_dir=paths[category]['autotrace'])

        # when process complete, delete input file
        logging.debug("Deleting %s" % s)
        os.remove(s)
    else:
        logging.error("Neuron is from unrecognized specimen category (%s)" % category)
    
def upload_all():
    parser = argparse.ArgumentParser()
    parser.add_argument('--nprocs', type=int, default=4)
    parser.add_argument('--log-level', default=logging.INFO)
    parser.add_argument('upload_dir')
    parser.add_argument('storage_dir')
    args = parser.parse_args()

    logging.basicConfig(level=args.log_level,
                        format='%(asctime)s %(process)d %(levelname)s: %(message)s', 
                        datefmt="%Y-%m-%d %H:%M:%S")


    pool = mp.Pool(processes=args.nprocs)
    
    paths = swank.storage_paths(args.storage_dir)

    for category, path_set in paths.items():
        for name, path in path_set.items():
            if not os.path.exists(path):
                os.makedirs(path)

    # perform upright on cells and copy to appropriate target project,
    #   then delete from upload directory
    swc_pattern = os.path.join(args.upload_dir, "*.swc")
    swcs = glob.glob(swc_pattern)

    logging.info("Found %d files in upload directory" % len(swcs))
    results = []
    for s in swcs:
        result = pool.apply_async(upload_swc, args=(s, paths))
        results.append(result)

    pool.close()

    for result in results:
        try:
            result.get()
        except InvalidMorphology as e:
            logging.error("Error uploading SWC: %s", e.message + " " + str(e.validation_errors))
        except Exception as e:
            logging.error("Error uploading SWC: %s", e.message)

    pool.join()

    # process mouse data
    logging.info("processing mouse")
    mouse.process_mouse_data(upright_to_process_dir=paths['Mouse']['upright_to_process'], 
                             upright_dir=paths['Mouse']['upright'],
                             non_upright_dir=paths['Mouse']['non_upright'], 
                             non_upright_to_process_dir=paths['Mouse']['non_upright_to_process'],
                             auto_dir=paths['Mouse']['autotrace'], 
                             features_dir=paths['Mouse']['features'], 
                             html_dir=paths['Mouse']['html'],
                             nprocs=args.nprocs)

    # TODO process human data
    logging.info("processing human")
    mouse.process_mouse_data(upright_to_process_dir=paths['Human']['upright_to_process'], 
                             upright_dir=paths['Human']['upright'],
                             non_upright_dir=paths['Human']['non_upright'], 
                             non_upright_to_process_dir=paths['Human']['non_upright_to_process'],
                             auto_dir=paths['Human']['autotrace'], 
                             features_dir=paths['Human']['features'], 
                             html_dir=paths['Human']['html'],
                             nprocs=args.nprocs)
    # TODO process lgn data

if __name__ == "__main__": upload_all()
