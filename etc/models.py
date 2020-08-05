import os

import numpy as np
from astropy import units as u
from synphot import units, SourceSpectrum, SpectralElement, specio
from synphot.spectrum import BaseUnitlessSpectrum, Empirical1D

__all__ = ['Site', 'Telescope', 'Instrument']

class Site:
    """Model for a site location and the atmosphere above it"""

    def __init__(self, name=None, altitude=None, latitude= None, longitude=None, **kwargs):
        self.name = name if name is not None else "Undefined"
        self.altitude = altitude * u.m if altitude is not None else altitude
        self.latitude = latitude
        self.longitude = longitude
        if 'transmission' in kwargs:
            modelclass = Empirical1D
            try:
                transmission = float(kwargs['transmission'])
                wavelengths = np.arange(300, 1501, 200) * u.nm
                throughput = len(wavelengths) * [transmission,]
                header = {}
            except ValueError:
                sky_file = os.path.expandvars(kwargs['transmission'])
                header, wavelengths, throughput = specio.read_spec(sky_file, wave_col='lam', flux_col='trans', wave_unit=u.nm,flux_unit=u.dimensionless_unscaled)
            self.transmission = BaseUnitlessSpectrum(modelclass, points=wavelengths, lookup_table=throughput, keep_neg=True, meta={'header': header})

    def __mul__(self, other):
        if isinstance(other, Telescope):
            other = other.reflectivity
        if isinstance(other, Instrument):
            other = other.transmission
        newcls = self
        newcls.transmission = self.transmission.__mul__(other)
        return newcls

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        newcls = self
        newcls.transmission = self.transmission.__truediv__(other)
        return newcls

    def __repr__(self):
        return "[{}({})]".format(self.__class__.__name__, self.name)

    def __str__(self):
        return "{}: lon={}, lat={}, altitude={}".format(self.name, self.longitude, self.latitude, self.altitude)


class Telescope:
    def __init__(self, name=None, size=0, area=0, num_mirrors=2, **kwargs):
        self.name = name if name is not None else "Undefined"
        self.size = size * u.m
        self.area_unit = u.m * u.m
        self.area = area * self.area_unit
        self.num_mirrors = num_mirrors

        modelclass = Empirical1D
        reflectivity = kwargs.get('reflectivity',  0.91)    # Default value based on average of bare Al over 300-1200nm
        try:
            reflectivity = float(reflectivity)
            wavelengths = np.arange(300, 1501, 200) * u.nm
            refl = len(wavelengths) * [reflectivity,]
            header = {}
        except ValueError:
            mirror_file = os.path.expandvars(kwargs['reflectivity'])
            print(mirror_file)
            header, wavelengths, refl = specio.read_ascii_spec(mirror_file, wave_unit=u.nm, flux_unit='%')

        mirror_se = BaseUnitlessSpectrum(modelclass, points=wavelengths, lookup_table=refl, keep_neg=True, meta={'header': header})
        # Assume all mirrors are the same reflectivity and multiply together
        self.reflectivity = mirror_se
        for x in range(0, self.num_mirrors-1):
            self.reflectivity *= mirror_se

    def tpeak(self, wavelengths=None):
        """Calculate :ref:`peak bandpass throughput <synphot-formula-tpeak>`.

        Parameters
        ----------
        wavelengths : array-like, `~astropy.units.quantity.Quantity`, or `None`
            Wavelength values for sampling.
            If not a Quantity, assumed to be in Angstrom.
            If `None`, ``self.waveset`` is used.

        Returns
        -------
        tpeak : `~astropy.units.quantity.Quantity`
            Peak bandpass throughput.

        """
        x = self.reflectivity._validate_wavelengths(wavelengths)
        return self.reflectivity(x).max()

    def __mul__(self, other):
        newcls = self
        newcls.reflectivity = self.reflectivity.__mul__(other)
        return newcls

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        newcls = self
        newcls.reflectivity = self.reflectivity.__truediv__(other)
        return newcls

    def __repr__(self):
        return "[{}({})]".format(self.__class__.__name__, self.name)

    def __str__(self):
        return "{} (M1: {} diameter, {} area; {} mirrors)".format(self.name, self.size.to(u.m), self.area.to(self.area_unit), self.num_mirrors)


class Instrument:
    def __init__(self, name=None, inst_type="IMAGER", **kwargs):
        _ins_types = ["IMAGER", "SPECTROGRAPH"]
        self.name = name if name is not None else "Undefined"
        self.inst_type = inst_type.upper() if inst_type.upper() in _ins_types else "IMAGER"

        # Defaults assume a "standard" imager with a single lens, 2 AR coatings
        # on the front and back surfaces and no mirrors
        self.num_ar_coatings = kwargs.get('num_ar_coatings', 2)
        self.num_lenses = kwargs.get('num_inst_lenses', 1)
        self.num_mirrors = kwargs.get('num_inst_mirrors', 0)

        # Fused silica (for the prism) and fused quartz (for the CCD window)
        # turn out to have the same transmission...
        self.lens_trans = kwargs.get('inst_lens_trans', 0.9)

        self.mirror_refl = kwargs.get('inst_mirror_refl', 0.9925)
        # Transmission/Reflection values of optical elements coating
        self.ar_coating = kwargs.get('inst_ar_coating_refl', 0.99)

        transmission = self._compute_transmission()
        wavelengths = np.arange(300, 1501, 200) * u.nm
        trans = len(wavelengths) * [transmission,]
        header = {}
        self.transmission = SpectralElement(Empirical1D, points=wavelengths, lookup_table=trans, keep_neg=True, meta={'header': header})

    def _compute_transmission(self):
        """This calculates the optical transmission of the instrument from lenses,
        mirrors and AR coatings. Assumes no/little wavelength dependence which
        is true for typical fused silica or quartz over most of optical/NIR regime
        see e.g.https://www.newport.com/n/optical-materials"""


        # Air-glass interfaces:
        throughput = self.ar_coating**self.num_ar_coatings
        # Transmissive optical elements
        throughput *= self.lens_trans**self.num_lenses
        # Reflective optical elements (Mirrors):
        throughput *= self.mirror_refl**self.num_mirrors

        return throughput



    def __repr__(self):
        return "[{}({})]".format(self.__class__.__name__, self.name)

