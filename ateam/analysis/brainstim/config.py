import os,sys
import json
import shutil


def set_dir_recursive(path,manifest):  
    '''
    Look recursively through the manifest and build the full paths
    '''

    path_split = path.split("/")
    dir_recursive = path_split[0]
    dir_static = "/".join(path_split[1:])
    if dir_recursive!='': # if not a BASE_DIR
        path = manifest[dir_recursive]
        path_recursive = set_dir_recursive(path,manifest)
        fullpath = "/".join([path_recursive,dir_static])
    else:
        fullpath=path
        
    return fullpath




def build(config_file):
    '''
    
    Set config
    
    '''
    
    with open(config_file, 'r') as f:
        conf=json.load(f)

    config_file_basename = os.path.basename(config_file)
    conf["config_file_name"] = config_file_basename
    conf["config_file_fullpath"] = os.path.abspath(config_file)

    manifest = build_manifest(conf)    
    recursive_insert(conf,manifest)
    
    return conf



def build_manifest(conf):

    '''
    Recursively look through the manifest to build the path variables (e.g.: $MY_DIR)
    '''
    manifest = conf["manifest"]
    
    for k,v in manifest.items():   # set paths in the manifest
        manifest[k] = set_dir_recursive(v,manifest)

        
    return manifest


def recursive_insert(d,manifest):

    '''
    Loop through the config and substitute the path variables (e.g.: $MY_DIR) with the values 
    from the manifest    
        
    '''
    for k,v in d.items():
        if isinstance(v,dict):
            d[k]=recursive_insert(v,manifest)
            
        else:
            if isinstance(v, basestring):

                path_split = v.split("/")
                dir_key = path_split[0]

                if dir_key!='' and dir_key[0]=='$':
                    path_split[0] = manifest[dir_key]
                
                d[k] = "/".join(path_split)
            
    return d


    
def copy(conf): 
    
    '''
    Copy config to output directory
    '''
    
    fullpath = conf['config_file_fullpath']
    print fullpath
    shutil.copy(fullpath, conf["manifest"]["$OUTPUT_DIR"])


def print_resolved(conf):

    '''
    Print the resolved config file to output directory
    '''
    name = conf['config_file_name'].split('.')
    name[-2] = name[-2] + "_resolved" # place just before extension
    file_name = '.'.join(name)
    fullpath = "/".join([conf["manifest"]["$OUTPUT_DIR"], file_name])

    with open(fullpath, 'w') as fp:
        json.dump(conf, fp, indent=4, separators=(',', ': '))  # print pretty
