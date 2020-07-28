
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
        self.components = []
