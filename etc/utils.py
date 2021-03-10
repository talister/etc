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
    header['source'] = source
    header['filename'] = filename
    element = BaseUnitlessSpectrum(Empirical1D, points=wavelengths, lookup_table=throughput, keep_neg=False, meta={'header': header})

    return element

def sptype_to_pickles_standard(sp_type):
    """Maps the passed <sp_type> e.g. 'F8V' to a Pickles standard star filename.
    None is returned if there is no match.
    References: http://www.stsci.edu/hst/observatory/crds/pickles_atlas.html
    Pickles (1998) (PASP 110, 863)"""

    # These are the main sequence dwarf, solar metallicity standards.
    mapping = { 'O5V' : 'pickles_1.fits',
                'O9V' : 'pickles_2.fits',
                'B0V' : 'pickles_3.fits',
                'B1V' : 'pickles_4.fits',
                'B3V' : 'pickles_5.fits',
                'B5V' : 'pickles_6.fits',
                'B6V' : 'pickles_6.fits',
                'B7V' : 'pickles_6.fits',
                'B8V' : 'pickles_7.fits',
                'B9V' : 'pickles_8.fits',
                'A0V' : 'pickles_9.fits',
                'A2V' : 'pickles_10.fits',
                'A3V' : 'pickles_11.fits',
                'A5V' : 'pickles_12.fits',
                'A7V' : 'pickles_13.fits',
                'F0V' : 'pickles_14.fits',
                'F2V' : 'pickles_15.fits',
                'F5V' : 'pickles_16.fits',
                'F6V' : 'pickles_18.fits',
                'F8V' : 'pickles_20.fits',
                'G0V' : 'pickles_23.fits',
                'G2V' : 'pickles_26.fits',
                'G5V' : 'pickles_27.fits',
                'G8V' : 'pickles_30.fits',
                'K0V' : 'pickles_31.fits',
                'K2V' : 'pickles_33.fits',
                'K3V' : 'pickles_34.fits',
                'K4V' : 'pickles_35.fits',
                'K5V' : 'pickles_36.fits',
                'K7V' : 'pickles_37.fits',
                'M0V' : 'pickles_38.fits',
                'M1V' : 'pickles_39.fits',
                'M2V' : 'pickles_40.fits',
                'M3V' : 'pickles_42.fits',
                'M4V' : 'pickles_43.fits',
                'M5V' : 'pickles_44.fits',
                'M6V' : 'pickles_45.fits',
                # Subgiants
                'B2IV' : 'pickles_46.fits'
            }

    return mapping.get(sp_type.upper(), None)
