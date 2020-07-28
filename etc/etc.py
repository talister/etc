
import toml
try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

from . import data
from .models import Site, Telescope, Instrument


class ETC(object):

    def __init__(self, components=None, config_file=None):
        PRESET_MODELS = toml.loads(pkg_resources.read_text(data, "FTN_FLOYDS.toml"))
        self.components = components if components is not None else []
        if config_file is None and len(self.components) == 0:
            component_config = PRESET_MODELS
        elif config_file is not None:
            component_config = toml.load(config_file)
        else:
            component_config = {}
        if len(component_config) > 0:
            self.create_components_from_config(component_config)

    def create_components_from_config(self, config):
        class_mapping = {'site' : Site, 'telescope' : Telescope, 'instrument' : Instrument }
        for component, component_config in config.items():
            class_name = class_mapping.get(component.lower(), None)
            if class_name:
#                print(class_name, component_config)
                output_component = class_name(**component_config)
                self.components.append(output_component)

    def __repr__(self):
        prefixstr = '<' + self.__class__.__name__ + ' '
        compstr = "->".join([repr(x) for x in self.components])
        return f'{prefixstr}{compstr}>'
