import os

import numpy as np
from astropy import units as u
from synphot import units, SourceSpectrum, SpectralElement, specio
from synphot.spectrum import BaseUnitlessSpectrum, Empirical1D

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
            self.transmission = SpectralElement(modelclass, points=wavelengths, lookup_table=throughput, keep_neg=True, meta={'header': header})

    def __mul__(self, other):

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

        self.reflectivity = SpectralElement(modelclass, points=wavelengths, lookup_table=refl, keep_neg=True, meta={'header': header})

    def __repr__(self):
        return "[{}({})]".format(self.__class__.__name__, self.name)

    def __str__(self):
        return "{} (M1: {} diameter, {} area; {} mirrors)".format(self.name, self.size.to(u.m), self.area.to(self.area_unit), self.num_mirrors)


class Instrument:
    def __init__(self, name=None, **kwargs):
        self.name = name if name is not None else "Undefined"

    def __repr__(self):
        return "[{}({})]".format(self.__class__.__name__, self.name)

