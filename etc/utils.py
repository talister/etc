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
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs import FITSFixedWarning
from astropy.table import QTable
from astropy.utils.exceptions import AstropyUserWarning
import numpy as np
from synphot import units, SourceSpectrum, SpectralElement, specio
from synphot.spectrum import BaseUnitlessSpectrum, Empirical1D

from .config import conf


class ETCError(Exception):
    pass


def read_element(filtername_or_filename, element_type='element', wave_units=u.nm, flux_units=units.THROUGHPUT):
    """Generic reader for optical elements, filters and others.
    The passed <filtername_or_filename> is first looked up in the `Conf` mapping
    dictionary for a match; if there is no match, it is assumed to be a filename
    or location specifier. These can be of 3 types:
    1. LCO Imaging Lab scanned optical element CSV files (contain 'LCO_' and '.csv' in the filename)
    2. SVO filter service references (contain 'http://svo')
    3. Local files (either ASCII or FITS)

    Checking on the read wavelengths is performed for local files:
    * if the first wavelength value is <100nm and the user didn't override the
    units through [wave_units], the wavelengths are assumed to be in, and are
    converted to, microns.
    * if the first wavelength value is >3000nm and the user didn't override the
    units through [wave_units], the wavelengths are assumed to be in, and are
    converted to, angstroms.
    """

    element_type = element_type.lower()
    filename =  conf.mapping.get(filtername_or_filename, None)
    if filename is None:
        filename = filtername_or_filename
    else:
        filename = filename()
        element_type = 'spectral_element'
    if 'LCO_' in filename.upper() and '.csv' in filename.lower():
        file_path = pkg_resources.files('etc.data').joinpath(os.path.expandvars(filename))
        source = "LCO iLab format"
        header, wavelengths, throughput = read_lco_filter_csv(file_path)
    elif 'http://svo' in filename.lower():
        source = "SVO filter service"
        header, wavelengths, throughput = specio.read_remote_spec(filename, wave_unit=u.AA, flux_unit=units.THROUGHPUT)
    else:
        source = "local file"
        file_path = os.path.expandvars(filename)
        if not os.path.exists(file_path):
            file_path = str(pkg_resources.files('etc.data').joinpath(filename))
        warnings.simplefilter('ignore', category = AstropyUserWarning)
        if filename.lower().endswith('fits') or filename.lower().endswith('fit'):
            wave_col = 'lam'
            flux_col = 'trans'
            if element_type == 'radiance':
                flux_col = 'flux'
            try:
                header, wavelengths, throughput = specio.read_spec(file_path, wave_col=wave_col, flux_col=flux_col, wave_unit=u.nm, flux_unit=flux_units)
            except KeyError:
                    # ESO-SM01 format; different column name for transmission and micron vs nm
                header, wavelengths, throughput = specio.read_spec(file_path, wave_col='lam', flux_col='flux', wave_unit=u.micron,flux_unit=flux_units)
        else:
            header, wavelengths, throughput = specio.read_ascii_spec(file_path, wave_unit=wave_units, flux_unit=flux_units)
        if wavelengths[0].value < 100.0 and wave_units == u.nm:
            # Small values seen, Convert to microns
            wavelengths = wavelengths.value * u.micron
        elif wavelengths[0].value > 3000.0 and wave_units == u.nm:
            # Large values seen, Convert to angstroms
            wavelengths = wavelengths.value * u.AA
        if element_type != 'spectrum' and element_type != 'radiance' and throughput.mean() > 1.5:
            # Test for mean throughput is above 1 to catch case where a throughput
            # fudge may be in the range ~1 to a few e.g. ESO Omegacam optics fudge
            # which goes to 3.2 and averages out to ~1.4
            throughput /= 100.0
            header['notes'] = 'Divided by 100.0 to convert from percentage'
    header['source'] = source
    header['filename'] = filename
    if element_type == 'spectrum' or element_type == 'radiance':
        # SourceSpectrum can't use the default units.THROUGHPUT so we need to
        # change to an assumed units.PHOTLAM (u.photon / (u.cm**2 * u.s * u.AA))
        # if nothing was passed by the user
        if flux_units == units.THROUGHPUT:
            if element_type == 'radiance':
                # Default units for ESO skycalc output (minus the arcsec^2)
                throughput = throughput.value * u.photon/u.s/u.m**2/u.um
            else:
                throughput = throughput.value * units.PHOTLAM
        element = SourceSpectrum(Empirical1D, points=wavelengths, lookup_table=throughput, keep_neg=False, meta={'header': header})
    elif element_type == 'spectral_element':
        element = SpectralElement(Empirical1D, points=wavelengths, lookup_table=throughput, keep_neg=False, meta={'header': header})
    else:
        element = BaseUnitlessSpectrum(Empirical1D, points=wavelengths, lookup_table=throughput, keep_neg=False, meta={'header': header})

    return element

def read_lco_filter_csv(csv_filter):
    """Reads filter transmission files in LCO Imaging Lab v1 format (CSV
    file with header and data)
    Returns an empty header dictionary and the wavelength and transmission columns"""

    table = QTable.read(csv_filter, format='ascii.csv', header_start=0, data_start=64)
    table.rename_column('ILDIALCT', 'Wavelength')
    table['Wavelength'].unit = u.nm
    table.rename_column('ilab_v1', 'Trans_measured')
    table['Trans_measured'].unit = u.dimensionless_unscaled
    table.rename_column('FITS/CSV file dialect', 'Trans_filtered')

    return {}, table['Wavelength'], table['Trans_measured']

def get_x_units(x_data):
    """finds wavelength units from x_data
       inputs: <xdata>: unitless wavelength data
       outputs: wavelength in Angstroms
    """
    x_mean = np.mean(x_data)

    # assuming visible to NIR range (~3000-10000A)
    if x_mean > 1000:
        x_units = u.AA  # (Angstroms)
    elif 100 < x_mean < 1000:
        x_units = u.nm
    elif .1 < x_mean < 10:
        x_units = u.micron
    else:
        print("Could not parse wavelength units from file. Assuming Angstroms")
        x_units = u.AA
    xxx = np.array(x_data)
    wavelength = (xxx * x_units).to(u.AA)
    return wavelength

def get_y_units(y_data, filename, header):
    """finds flux/reflectance units
       inputs: y_data, spectrum file
       outputs: scaled flux with Units
    """
    y_factor = 1
    if 'BUNIT' in header:
        try:
            y_units = u.Unit(header['BUNIT'])
        except ValueError:
            print("Could not parse flux units from header.")

    elif "ctiostan" in filename and '.dat' in filename:  # from ESO aaareadme.ctio
        y_units = u.erg/(u.cm**2)/u.s/u.AA
        y_factor = 10**16

    elif .001 < np.median(y_data) < 10:  # Probably Normalized
        y_units = u.def_unit("Normalized_Reflectance", u.dimensionless_unscaled)

    elif '_2df_ex.fits' in filename:  # from FLOYDS
        y_factor = 10**20
        y_units = u.erg/(u.cm**2)/u.s/u.AA

    else:
        print("Could not parse flux units from file. Assuming erg/cm^2/s/A")
        y_units = u.erg/(u.cm**2)/u.s/u.AA

    yyy = np.array(y_data)
    flux = ((yyy / y_factor) * y_units)
    return flux

def read_eso_spectra(filepath):

    try:
        hdul = fits.open(filepath)
    except FileNotFoundError:
        print("Cannot find file {}".format(filepath))
        return None

    data = hdul[0].data
    hdr = hdul[0].header

    yyy = data
    if hdr.get('NAXIS3') == 4:
        # FLOYDS merged data, slice/bandid 1 has the extracted spectrum
        yyy = data[0][0]
    warnings.simplefilter('ignore', category=FITSFixedWarning)
    w = WCS(hdr, naxis=1, relax=False, fix=False)
    lam = w.wcs_pix2world(np.arange(len(yyy)), 0)[0]

    wavelength = get_x_units(lam)
    flux = get_y_units(yyy, filepath, hdr)

    source_spec = SourceSpectrum(Empirical1D, points=wavelength, lookup_table=flux,
                   keep_neg=True, meta={'header': hdr})
    return source_spec

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

def percentage_difference(v1, v2):
    """Calculates the percentage difference between value <v1> and <v2>
    The result is returned as an astropy percentage Quantity"""

    return np.abs(v1-v2)/((v1+v2)/2.)*100*u.percent
