import os
import sys
import warnings

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
from astropy.utils.exceptions import AstropyUserWarning
from synphot import units, SourceSpectrum, SpectralElement, specio
from synphot.spectrum import BaseUnitlessSpectrum, Empirical1D

from .config import conf


class ETCError(Exception):
    pass


def read_element(filtername_or_filename, wave_units=u.nm, flux_units=units.THROUGHPUT):
    filename =  conf.mapping.get(filtername_or_filename, None)
    if filename is None:
        filename = filtername_or_filename
    if 'LCO_' in filename.upper() and '.csv' in filename.lower():
        file_path = pkg_resources.files('etc.data').joinpath(os.path.expandvars(filename))
        source = "LCO iLab format"
        header, wavelengths, throughput  = self._read_lco_filter_csv(file_path)
    elif 'http://svo' in filename.lower():
        source = "SVO filter service"
        header, wavelengths, throughput = specio.read_remote_spec(filename, wave_unit=u.AA, flux_unit=units.THROUGHPUT)
    else:
        source = "local file"
        file_path = str(pkg_resources.files('etc.data').joinpath(os.path.expandvars(filename)))
        warnings.simplefilter('ignore', category = AstropyUserWarning)
        if filename.lower().endswith('fits') or filename.lower().endswith('fit'):
            try:
                header, wavelengths, throughput = specio.read_spec(file_path, wave_col='lam', flux_col='trans', wave_unit=u.nm, flux_unit=u.dimensionless_unscaled)
            except KeyError:
                    # ESO-SM01 format; different column name for transmission and micron vs nm
                header, wavelengths, throughput = specio.read_spec(file_path, wave_col='lam', flux_col='flux', wave_unit=u.micron,flux_unit=u.dimensionless_unscaled)
        else:
            header, wavelengths, throughput = specio.read_ascii_spec(file_path, wave_unit=wave_units, flux_unit=flux_units)
        if throughput.mean() > 1.0:
            throughput /= 100.0
            header['notes'] = 'Divided by 100.0 to convert from percentage'

    element = BaseUnitlessSpectrum(Empirical1D, points=wavelengths, lookup_table=throughput, keep_neg=False, meta={'header': header})

    return element
