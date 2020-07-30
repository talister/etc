import os

import numpy as np
from astropy import units as u
from synphot import units, SourceSpectrum, SpectralElement, specio
from synphot.spectrum import BaseUnitlessSpectrum, Empirical1D, ConstFlux1D

class Site(SpectralElement):
    """Model for a site location and the atmosphere above it"""

    def __init__(self, name=None, altitude=None, latitude= None, longitude=None, *args, **kwargs):
        self.name = name if name is not None else "Undefined"
        self.altitude = altitude * u.m if altitude is not None else altitude
        self.latitude = latitude
        self.longitude = longitude
        self.meta = {}
        if 'transmission' in kwargs:
            modelclass = Empirical1D
            try:
                transmission = float(kwargs['transmission'])
                print("Is float")
                wavelengths = np.arange(300, 1501, 200) * u.nm
                throughput = len(wavelengths) * [transmission,]
                header = {}
            except ValueError:
                sky_file = os.path.expandvars(kwargs['transmission'])
                print("Is filename", sky_file)
                header, wavelengths, throughput = specio.read_spec(sky_file, wave_col='lam', flux_col='trans', wave_unit=u.nm,flux_unit=u.dimensionless_unscaled)
            super().__init__(modelclass, points=wavelengths, lookup_table=throughput, keep_neg=True, meta={'header': header})
            print("model=",self._model)
        # else:
            # wavelengths, lookup_table = self._get_arrays(self.waveset)
            # print("model=",self._model)
            # super().__init__(self._model, points=wavelengths, lookup_table=lookup_table, *args, **kwargs)

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

    def __repr__(self):
        return "[{}({})]".format(self.__class__.__name__, self.name)

    def __str__(self):
        return "{} (M1: {} diameter, {} area; {} mirrors)".format(self.name, self.size.to(u.m), self.area.to(self.area_unit), self.num_mirrors)


class Instrument:
    def __init__(self, name=None, **kwargs):
        self.name = name if name is not None else "Undefined"

    def __repr__(self):
        return "[{}({})]".format(self.__class__.__name__, self.name)

