import os
import sys
import numbers
import warnings
from collections import OrderedDict

try:
    if sys.version_info >= (3, 9):
        # This exists in 3.8 but a different API
        import importlib.resources as pkg_resources
    else:
        raise ImportError
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

import numpy as np
from astropy import units as u
from astropy.table import QTable
from astropy.io.ascii import InconsistentTableError
from astropy.utils.exceptions import AstropyUserWarning
from synphot import units, SourceSpectrum, SpectralElement, specio
from synphot.spectrum import BaseUnitlessSpectrum, Empirical1D
from synphot.observation import Observation

from .config import conf
from .utils import ETCError, read_element

__all__ = ['Site', 'Telescope', 'Instrument']


class Site:
    """Model for a site location and the atmosphere above it"""

    def __init__(self, name=None, altitude=None, latitude= None, longitude=None, **kwargs):
        self.name = name if name is not None else "Undefined"
        self.altitude = altitude * u.m if altitude is not None else altitude
        self.latitude = latitude * u.deg if latitude is not None else latitude
        self.longitude = longitude * u.deg if longitude is not None else longitude
        if 'transmission' in kwargs:
            modelclass = Empirical1D
            try:
                transmission = float(kwargs['transmission'])
                wavelengths = np.arange(300, 1501, 1) * u.nm
                throughput = len(wavelengths) * [transmission,]
                header = {}
            except ValueError:
                # XXX Replace with read_element()?
                sky_file = str(pkg_resources.files('etc.data').joinpath(os.path.expandvars(kwargs['transmission'])))
                try:
                    header, wavelengths, throughput = specio.read_spec(sky_file, wave_col='lam', flux_col='trans', wave_unit=u.nm,flux_unit=u.dimensionless_unscaled)
                except KeyError:
                    # ESO-SM01 format; different column name for transmission and micron vs nm
                    header, wavelengths, throughput = specio.read_spec(sky_file, wave_col='lam', flux_col='flux', wave_unit=u.micron,flux_unit=u.dimensionless_unscaled)
                except InconsistentTableError:
                    # ASCII format ?
                    header, wavelengths, throughput = specio.read_spec(sky_file, wave_unit=u.nm,flux_unit=u.dimensionless_unscaled)
            self.transmission = BaseUnitlessSpectrum(modelclass, points=wavelengths, lookup_table=throughput, keep_neg=False, meta={'header': header})
        if 'sky_mag' in kwargs:
            self.sky_mags = kwargs['sky_mag']
        else:
            file_path = pkg_resources.files('etc.data').joinpath(os.path.expandvars(conf.sky_brightness_file))
            self.sky_mags_table = self._read_skybrightness_file(file_path)
            self.sky_mags = []
            if self.sky_mags_table:
                self.sky_mags = dict(zip(self.sky_mags_table.colnames, self.sky_mags_table[0]))

    def sky_spectrum(self, filtername='V'):
        waveset = np.arange(3000,12001,1) * u.AA
        sky_flux = np.empty(len(waveset))
        sky_flux.fill(self._photon_rate(filtername))
        sky_flux = u.Quantity(sky_flux, unit=units.PHOTLAM)
        sky = SourceSpectrum(Empirical1D, points=waveset, lookup_table=sky_flux)

        return sky

    def rebin_transmission(self, wavenew):
        """Rebins a potentially higher resolution atmospheric transmission
        spectrum to a new wavelength grid <wavenew>.
        Returns a new BaseUnitlessSpectrum constructed from the new passed
        wavelength grid and the resampled transmission"""

        wave = self.transmission.waveset
        specin = self.transmission(wave)
        # Have to pretend transmission is in a flux unit to make a SourceSpectrum
        # so drop here and put back at the end
        orig_unit = specin.unit
        spec = SourceSpectrum(Empirical1D, points=wave, lookup_table=specin.value)
        f = np.ones(len(wave))
        filt = SpectralElement(Empirical1D, points=wave, lookup_table=f)
        obs = Observation(spec, filt, binset=wavenew, force='taper')

        new_trans = SpectralElement(Empirical1D, points=wavenew, lookup_table=obs.binflux.value * orig_unit)

        return new_trans

    def tpeak(self, wavelengths=None):
        """Calculate peak atmospheric transmission`.

        Parameters
        ----------
        wavelengths : array-like, `~astropy.units.quantity.Quantity`, or `None`
            Wavelength values for sampling.
            If not a Quantity, assumed to be in Angstrom.
            If `None`, ``self.waveset`` is used.

        Returns
        -------
        tpeak : `~astropy.units.quantity.Quantity`
            Peak atmospheric transmission.

        """
        x = self.transmission._validate_wavelengths(wavelengths)
        return self.transmission(x).max()

    def _extinction_to_transmission(self, extinction, airmass=1.0):
        """Transform a passed extinction coefficient into a fractional transmission"""

        transmission = 10**(-0.4*extinction * airmass)

        return transmission

    def _transmission_to_extinction(self, transmission, airmass=1.0):
        """Transform a passed fractional transmission into an extinction coefficient"""

        extinct = np.log10(transmission) / (-0.4 * airmass)

        return extinct

    def _photon_rate(self, filtername='V'):
        """Calculates the photon rate in photons/s/cm^2/AA for mag=0 for the
        passed [filtername] (defaults to V)
        Values originally from SIGNAL and hence Bessell (1979 PASP 91, 589), now
        updated to Bessell, Castelli & Plez (1998 A&A 333, 231) and Fukugita et al. 1996
        for SDSS/PanSTARRS
        """

        flux_janskys = {'U': 1790, 'B': 4063, 'V' : 3636, 'R' : 3064, 'I' : 2416, 'Z' : 2200,
                    'u' : 3631, 'g': 3631, 'r': 3631, 'i': 3631, 'z': 3631,
                    'up' : 3631, 'gp': 3631, 'rp': 3631, 'ip': 3631, 'zp': 3631, 'w' : 3631}
        flux_mag0_Jy = flux_janskys[filtername] * u.Jy
        wavelength = self._map_filter_to_wavelength(filtername)
        m_0 = flux_mag0_Jy.to(u.photon / u.cm**2 / u.s / u.angstrom, equivalencies=u.spectral_density(wavelength))

        return m_0

    def _map_filter_to_wavelength(self, filtername='V'):
        """Maps the given [filtername] (defaults to 'V' for Bessell-V') to a wavelength
        which is returned as an AstroPy Quantity in angstroms"""

        filter_cwave = {'U': 3600, 'B': 4300, 'V' : 5500, 'R' : 6500, 'I' : 8200, 'Z' : 9500,
                        'u' : 3675, 'g' : 4763, 'rp' : 6204, 'ip' : 7523, 'zp' : 8660, 'z' : 9724,
                        'gp' : 4810, 'rp' : 6170, 'ip' : 7520, 'zp' : 8660, 'w' : 6080}
        wavelength = filter_cwave[filtername] * u.angstrom

        return wavelength

    def _read_skybrightness_file(self, sky_file):
        table = None
        try:
            table = QTable.read(sky_file, format='ascii.commented_header')
        except FileNotFoundError:
            print("Cannot find", sky_file)
        except:
            raise

        return table

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
        mirrors = []
        # Default value based on average of coated Al over 300-1200nm and also matches SIGNAL
        reflectivity = kwargs.get('reflectivity',  0.85)
        try:
            reflectivity = float(reflectivity)
            wavelengths = np.arange(300, 1501, 1) * u.nm
            refl = len(wavelengths) * [reflectivity,]
            header = {}
            mirror_se = BaseUnitlessSpectrum(modelclass, points=wavelengths, lookup_table=refl, keep_neg=True, meta={'header': header})
            # Assume all mirrors are the same reflectivity and copy
            for x in range(self.num_mirrors):
                mirrors.append(mirror_se)
        except ValueError:
            # single filename
            component =  kwargs['reflectivity']
            mirror_se =read_element(component)
            for x in range(self.num_mirrors):
                mirrors.append(mirror_se)
        except TypeError:
            # List of filename components
            telescope_components = kwargs['reflectivity']
            for component in telescope_components:
                mirror_se = read_element(component)
                mirrors.append(mirror_se)
        # Multiply mirror reflectivities together
        self.reflectivity = mirrors[0]
        for x in range(0, self.num_mirrors-1):
            self.reflectivity *= mirrors[x+1]

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
    adc_error = np.sqrt(0.289) * (u.adu / u.pixel)

    def __init__(self, name=None, inst_type="IMAGER", **kwargs):
        _ins_types = ["IMAGER", "SPECTROGRAPH"]
        self.name = name if name is not None else "Undefined"
        self.inst_type = inst_type.upper() if inst_type.upper() in _ins_types else "IMAGER"

        # Defaults assume a "standard" imager with a single lens, 2 AR coatings
        # on the front and back surfaces and no mirrors
        self.num_ar_coatings = kwargs.get('num_ar_coatings', 2)
        self.num_lenses = kwargs.get('num_inst_lenses', 1)
        self.num_mirrors = kwargs.get('num_inst_mirrors', 0)

        # Fused silica (common lens material) and fused quartz (common for CCD windows)
        # turn out to have the same transmission...
        self.lens_trans = kwargs.get('inst_lens_trans', 0.93)

        self.mirror_refl = kwargs.get('inst_mirror_refl', 0.9925)
        # Transmission/Reflection values of optical elements coating
        self.ar_coating = kwargs.get('inst_ar_coating_refl', 0.99)

        trans_components = kwargs.get('trans_components',  None)
        wavelengths = np.arange(300, 1501, 1) * u.nm
        if trans_components:
            trans = len(wavelengths) * [1.0,]
            for comp_name in trans_components.split(","):
                print(comp_name)
                element = read_element(comp_name)
                trans *= element(wavelengths)
        else:
            print("Computing transmission from elements")
            transmission = self._compute_transmission()
            trans = len(wavelengths) * [transmission,]
        header = {}
        self.transmission = SpectralElement(Empirical1D, points=wavelengths, lookup_table=trans, keep_neg=True, meta={'header': header})

        self.filterlist = kwargs.get('filterlist', [])
        self.filterset = OrderedDict()
        for filtername in self.filterlist:
            if filtername not in self.filterset:
                self.filterset[filtername] = self.set_bandpass_from_filter(filtername)

        ccd_qe = kwargs.get('ccd_qe', 0.9)
        if not isinstance(ccd_qe, (u.Quantity, numbers.Number)):
            file_path = os.path.expandvars(ccd_qe)
            if not os.path.exists(file_path):
                file_path = str(pkg_resources.files('etc.data').joinpath(ccd_qe))
            header, wavelengths, throughput = specio.read_ascii_spec(file_path, wave_unit=u.nm, flux_unit=units.THROUGHPUT)
            if throughput.mean() > 1.0:
                throughput /= 100.0
                header['notes'] = 'Divided by 100.0 to convert from percentage'
            header['filename'] = ccd_qe
            self.ccd_qe = BaseUnitlessSpectrum(Empirical1D, points=wavelengths, lookup_table=throughput, keep_neg=False, meta={'header': header})
        else:
            self.ccd_qe = ccd_qe
        self.ccd_gain = kwargs.get('ccd_gain', 1) * (u.electron / u.adu)
        self.ccd_readnoise = kwargs.get('ccd_readnoise', 0) * (u.electron / u.pix)
        ccd_pixsize = kwargs.get('ccd_pixsize', 0)
        try:
            ccd_pixsize_units = u.Unit(kwargs.get('ccd_pixsize_units', 'micron'))
        except ValueError:
            ccd_pixsize_units = u.micron
        self.ccd_pixsize = (ccd_pixsize * ccd_pixsize_units).to(u.micron)
        self.ccd_xpixels = kwargs.get('ccd_xpixels', 0)
        self.ccd_ypixels = kwargs.get('ccd_ypixels', 0)

        fwhm = kwargs.get('fwhm', 1)
        try:
            fwhm_units = u.Unit(kwargs.get('fwhm_units', 'arcsec'))
        except ValueError:
            fwhm_units = u.arcsec
        try:
            # Already a Quantity
            self.fwhm = fwhm.to(u.arcsec)
        except AttributeError:
            self.fwhm = fwhm * fwhm_units

        focal_scale = kwargs.get('focal_scale', 0)
        try:
            focal_scale_units = u.Unit(kwargs.get('focal_scale_units', 'arcsec/mm'))
        except ValueError:
            focal_scale_units = u.arcsec/u.mm
        self.focal_scale = focal_scale * focal_scale_units
        if self.ccd_pixsize !=0 and self.focal_scale != 0:
            self.ccd_pixscale = self.focal_scale.to(u.arcsec/u.mm) * self.ccd_pixsize.to(u.mm)

    @property
    def is_imager(self):
        return True if self.inst_type == 'IMAGER' else False

    def ccd_fov(self, fov_units=u.arcsec):
        """Computes the CCD's field of view and returns a tuple of Quantity's
        Uses `self.ccd_pixsize`, `self.focal_scale` and `self.ccd_x/ypixels`
        to compute the size"""

        xsize=ysize=0*u.arcsec
        if self.ccd_pixsize is not None and self.focal_scale is not None and\
            self.ccd_xpixels is not None and self.ccd_ypixels is not None:
            xsize = self.ccd_xpixels * self.ccd_pixsize * self.focal_scale
            ysize = self.ccd_ypixels * self.ccd_pixsize * self.focal_scale
        return xsize.decompose().to(fov_units), ysize.decompose().to(fov_units)

    def _read_lco_filter_csv(self, csv_filter):
        """Reads filter transmission files in LCO Imaging Lab v1 format (CSV
        file with header and data)
        Returns an empty header dictionary and the wavelength and trensmission columns"""

        table = QTable.read(csv_filter, format='ascii.csv', header_start=0, data_start=64)
        table.rename_column('ILDIALCT', 'Wavelength')
        table['Wavelength'].unit = u.nm
        table.rename_column('ilab_v1', 'Trans_measured')
        table['Trans_measured'].unit = u.dimensionless_unscaled
        table.rename_column('FITS/CSV file dialect', 'Trans_filtered')

        return {}, table['Wavelength'], table['Trans_measured']


    def set_bandpass_from_filter(self, filtername):
        """Loads the specified <filtername> from the transmission profile file
        which is mapped via the etc.config.Conf() items.
        Returns a SpectralElement instance for the filter profile
        """

        if len(filtername) == 2 and filtername[1] == 'p':
            filtername = filtername[0]

        filename =  conf.mapping.get(filtername, None)
        if filename is None:
            raise ETCError('Filter name {0} is invalid.'.format(filtername))
        if 'LCO_' in filename().upper() and '.csv' in filename().lower():
            file_path = pkg_resources.files('etc.data').joinpath(os.path.expandvars(filename()))
            source = "LCO iLab format"
            header, wavelengths, throughput  = self._read_lco_filter_csv(file_path)
        elif 'http://svo' in filename().lower():
            source = "SVO filter service"
            header, wavelengths, throughput = specio.read_remote_spec(filename(), wave_unit=u.AA, flux_unit=units.THROUGHPUT)
        else:
            source = "local file"
            file_path = pkg_resources.files('etc.data').joinpath(os.path.expandvars(filename()))
            warnings.simplefilter('ignore', category = AstropyUserWarning)
            header, wavelengths, throughput = specio.read_ascii_spec(file_path, wave_unit=u.nm, flux_unit=units.THROUGHPUT)
            if throughput.mean() > 1.0:
                throughput /= 100.0
                header['notes'] = 'Divided by 100.0 to convert from percentage'
        print("Reading from {} for {}".format(source, filtername))

        header['filename'] = filename
        header['descrip'] = filename.description
        meta = {'header': header, 'expr': filtername}

        return SpectralElement(Empirical1D, points=wavelengths, lookup_table=throughput, meta=meta)

    def slit_vignette(self, slit_width=1*u.arcsec):
        """Compute the fraction of light entering the slit of width <slit_width>
        for an object described by a FWHM of <self.fwhm>
        In the case of imaging mode, 1.0 is always returned.

        The code is taken from the IAC/Chris Benn's SIGNAL code (`slitvign` routine):
        http://www.ing.iac.es/Astronomy/instruments/signal/help.html
        http://www.ing.iac.es/Astronomy/instruments/signal/signal_code.html
        which in turn is an approximation (<5% error) to the numerical simulation
        from the `light_in_slit` code (http://www.ing.iac.es/~crb/misc/lightinslit.f)
        Assuming a Gaussian/Normal distribution and ignoring differential refraction,
        this can also be calculated analytically as:
            Phi(n) - Phi(-n)
        where Phi(x) is the normal function cumulative distribution function (e.g.
        `scipy.stats.norm.cdf()`  and `n` is the slit width sigma/2 (since the
        distribution is symmetric about 0).
        """

        vign = 1.0

        if self.is_imager is False:
            # Spectroscopy
            try:
                ratio = slit_width.to(u.arcsec) / self.fwhm.to(u.arcsec)
            except AttributeError:
                ratio = slit_width / self.fwhm.value

            if ratio < 0.76:
                vign = 0.868*ratio
            if 0.76 <= ratio < 1.40:
                vign = 0.37+0.393*ratio
            if 1.40 <= ratio < 2.30:
                vign = 1.00-0.089*(2.3-ratio)
            if ratio >= 2.3:
                vign = 1.0

        return vign

    def throughput(self, filtername):
        """Returns the total throughput of optics+filter/grating+CCD"""

        if filtername not in self.filterlist or filtername not in self.filterset:
            raise ETCError('Filter name {0} is invalid.'.format(filtername))

        return self.filterset[filtername] * self.transmission * self.ccd_qe

    def _compute_transmission(self):
        """This calculates the optical transmission of the instrument from lenses,
        mirrors and AR coatings. Assumes no/little wavelength dependence which
        is true for typical fused silica or quartz over most of optical/NIR regime
        see e.g. https://www.newport.com/n/optical-materials"""

        # Air-glass interfaces:
        throughput = self.ar_coating**self.num_ar_coatings
        # Transmissive optical elements
        throughput *= self.lens_trans**self.num_lenses
        # Reflective optical elements (Mirrors):
        throughput *= self.mirror_refl**self.num_mirrors

        return throughput

    def __repr__(self):
        return "[{}({})]".format(self.__class__.__name__, self.name)

