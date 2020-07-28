
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
            self._create_components_from_config(component_config)

    def _create_components_from_config(self, config):
        class_mapping = {'site' : Site, 'telescope' : Telescope, 'instrument' : Instrument }
        for component, component_config in config.items():
            class_name = class_mapping.get(component.lower(), None)
            if class_name:
#                print(class_name, component_config)
                output_component = class_name(**component_config)
                self.components.append(output_component)

    @property
    def site(self):
        sites = [x for x in self.components if isinstance(x, Site)]
        return sites[0] if len(sites) == 1 else None

    @property
    def telescope(self):
        tels = [x for x in self.components if isinstance(x, Telescope)]
        return tels[0] if len(tels) == 1 else None

    @property
    def instrument(self):
        insts = [x for x in self.components if isinstance(x, Instrument)]
        return insts[0] if len(insts) == 1 else None

    def __repr__(self):
        prefixstr = '<' + self.__class__.__name__ + ' '
        compstr = "->".join([repr(x) for x in self.components])
        return f'{prefixstr}{compstr}>'
