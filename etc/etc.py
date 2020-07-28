
import toml
try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

from . import data


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
        for component, component_config in config.items():
            self.components.append(component)
