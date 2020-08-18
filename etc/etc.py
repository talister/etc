import sys
import toml
try:
    if sys.version_info >= (3, 9):
        # This exists in 3.8 but a different API
        import importlib.resources as pkg_resources
    else:
        raise ImportError
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

from astropy import units as u

from . import data
from .models import Site, Telescope, Instrument


class ETC(object):

    _internal_wave_unit = u.nm

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

        self.combined = None

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

    def _do_plot(self, waves, thru, title='', left=300*u.nm, right=1200*u.nm, bottom=0.0, top=None):
        """Plot worker.

        Parameters
        ----------
        waves, trans : `~astropy.units.quantity.Quantity`
            Wavelength and throughput to plot.

        kwargs
            See :func:`plot`.

        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print('No matplotlib installation found; plotting disabled as a result.')
            return

        fig, ax = plt.subplots()
        ax.plot(waves, thru)

        # Custom wavelength limits
        if left is not None:
            if isinstance(left, u.Quantity):
                left_val = left.to(self._internal_wave_unit).value
            else:
                left_val = left
            ax.set_xlim(left=left_val)
        if right is not None:
            if isinstance(right, u.Quantity):
                right_val = right.to(self._internal_wave_unit).value
            else:
                right_val = right
            ax.set_xlim(right=right_val)

        # Custom flux/throughput limit
        if bottom is not None:
            ax.set_ylim(bottom=bottom)
        if top is not None:
            ax.set_ylim(top=top)

        wave_unit = waves.unit
        ax.set_xlabel('Wavelength ({0})'.format(wave_unit))

        yu = thru.unit
        if yu is u.dimensionless_unscaled:
            ax.set_ylabel('Unitless')
        else:
            ax.set_ylabel('Flux ({0})'.format(yu))

        if title:
            ax.set_title(title)

        # Turn on minorticks on both axes
        ax.minorticks_on()

        plt.draw()

    def plot(self, **kwargs):
        self._create_combined()
        waves, trans = self.combined._get_arrays(None)
        self._do_plot(waves.to(self._internal_wave_unit), trans, **kwargs)

    def _create_combined(self):
        if self.combined is None:
            self.combined = self.site.transmission * self.telescope.reflectivity * self.instrument.transmission * self.instrument.ccd

    def __repr__(self):
        prefixstr = '<' + self.__class__.__name__ + ' '
        compstr = "->".join([repr(x) for x in self.components])
        return f'{prefixstr}{compstr}>'
