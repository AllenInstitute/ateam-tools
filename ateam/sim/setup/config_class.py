import json
import os.path
import collections

class ConfigBuilder(dict):

    def __init__(self, *args, **kwargs):
        super(ConfigBuilder, self).__init__(*args, **kwargs)
        self._path = ""

    def save_copy_as(self, output_path):
        with open(output_path, 'w') as fp:
            json.dump(self, fp, indent=2)

    def save(self):
        self.save_copy_as(self.path)

    def save_as(self, output_path):
        self._path = output_path
        self.save()
    
    @staticmethod
    def update_nesteddict(d, *dicts, **kwargs):
        for udict in dicts + (kwargs,):
            for k, v in udict.iteritems():
                if isinstance(v, collections.Mapping):
                    d[k] = ConfigBuilder.update_nesteddict(d.get(k, {}), v)
                else:
                    d[k] = v
        return d

    def update_nested(self, *dicts, **kwargs):
        return self.update_nesteddict(self, *dicts, **kwargs)

    def fix_components_path(self):
        components_dir = self.get('manifest', {}).get('$COMPONENTS_DIR')
        components_path = os.path.normpath(os.path.join(self.folder, components_dir))
        self.update_nested({'manifest':{'$COMPONENTS_DIR': components_path}})
        return components_path

    @property
    def path(self):
        return os.path.abspath(self._path)

    @property
    def folder(self):
        return os.path.dirname(self.path)

    @classmethod
    def from_json(cls, config_file):
        conf = cls(json.load(open(config_file, 'r')))
        conf._path = config_file
        return conf

    @classmethod
    def load_template(cls, config_template, output_path, shared_components=True):
        """ Copies config file with absolute paths, loading it to resolve paths. """
        # TODO: check rest of manifest etc?
        conf = cls.from_json(config_template)
        if shared_components:
            conf.fix_components_path()
        conf.save_as(output_path)
        return conf