import os
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

import numpy as np
from astropy import units as u
from synphot import SourceSpectrum, SpectralElement, units
from synphot.models import Empirical1D
from synphot.observation import Observation

from . import data
from .models import Site, Telescope, Instrument
from .config import conf


class ETC(object):

    _internal_wave_unit = u.nm
    _sun_raw = SourceSpectrum.from_file(os.path.expandvars(conf.sun_file))
    _vega = SourceSpectrum.from_file(os.path.expandvars(conf.vega_file))
    _V_band = SpectralElement.from_filter('johnson_v')

    def __init__(self, config_file=None, components=None):
        PRESET_MODELS = toml.loads(pkg_resources.read_text(data, "FTN_FLOYDS.toml"))
        self.components = components if components is not None else []
        if config_file is None and len(self.components) == 0:
            component_config = PRESET_MODELS
        elif config_file is not None:
            if type(config_file) == dict:
                component_config = config_file
            else:
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

    def photons_from_source(self, V_mag, filtername):
        sun = self._vega.normalize(V_mag * units.VEGAMAG, self._V_band, vegaspec=self._vega)

        source_spec = sun
        self._create_combined()
        waves, thru = self.combined._get_arrays(None)
        filter_waves, filter_trans = self.instrument.filterset[filtername]._get_arrays(waves)
        spec_elements = SpectralElement(Empirical1D, points=filter_waves, lookup_table = filter_trans * thru)

        # get the synphot observation object
        synphot_obs = Observation(source_spec, spec_elements, force='taper')

        countrate = synphot_obs.countrate(area=self.telescope.area)

        return countrate

    def ccd_snr(self, exp_time, V_mag, filtername, npix=1 * u.pixel,
                         n_background=np.inf * u.pixel,
                         background_rate=0 * (u.ph / u.pixel / u.s),
                         sky_mag=None,
                         darkcurrent_rate=0 * (u.ph / u.pixel / u.s)):
        """Calculate the SNR in a given exposure time <exp_time> for a given <V_mag>
        """

        try:
            exp_time = exp_time.to(u.s)
        except AttributeError:
            exp_time = exp_time * u.s

        # Get countrate per second for this magnitude in the filter
        countrate = self.photons_from_source(V_mag, filtername)

        # define countrate to be in photons/sec, which can be treated as the
        # same as (photo)electrons for CCDs but are not technically convertible in
        # astropy.units:
        countrate = countrate.value * (u.photon / u.s)

        print("Counts/s=", countrate)

        readnoise = self.instrument.ccd_readnoise.value * (u.photon/u.pixel)
        readnoise_sq = readnoise.value**2 * (u.photon/u.pixel)
        gain_err = self._get_shotnoise(self.instrument.ccd_gain * self.instrument.adc_error)
        print("Readnoise=", readnoise, readnoise_sq, gain_err)
        if sky_mag is not None:
            # Compute countrate from given sky magnitude.
            # XXX need to refactor photons_from_source() to do this also
            sky_filtername = filtername
            if '::' in filtername:
                sky_filtername = filtername.split('::')[1]
            sky = self.site.sky_spectrum(sky_filtername)
            print("  Original sky=" ,sky(sky.avgwave()))
            sky2 = sky.normalize(sky_mag*units.VEGAMAG, self._V_band, vegaspec=self._vega)
            print("Normalized sky=", sky2(sky2.avgwave()), sky(sky.avgwave())*10**(-sky_mag/2.5))
            self._create_combined()
            waves, thru = self.combined._get_arrays(None)
            obs_filter = self.instrument.filterset[filtername]
            filter_waves, filter_trans = obs_filter._get_arrays(waves)
            print('thruput=', thru.max(), self.telescope.area.to(units.AREA)*thru.max())
            spec_elements = SpectralElement(Empirical1D, points=filter_waves, lookup_table = thru * filter_trans)

            # get the synphot observation object
            sky_obs = Observation(sky2, spec_elements, force='taper')

            sky_countrate = sky_obs.countrate(area=self.telescope.area)
            # Actually a count rate per arcsec^2 as the original magnitude is
            # also per sq. arcsec. Also set to photons/s to match signal from object
            sky_countrate = sky_countrate.value * (u.photon/u.s/u.arcsec**2)
            print("Sky (photons/AA/arcsec^2)=", sky_countrate, sky_countrate/obs_filter.equivwidth())
            sky_per_pixel = sky_countrate * self.instrument.ccd_pixscale * self.instrument.ccd_pixscale
            sky_per_pixel /= u.pixel
            print("Sky (photons/pixel)=", sky_per_pixel)
            background_rate = sky_per_pixel

        fwhm = self.instrument.fwhm
        pixel_scale = self.instrument.ccd_pixscale
        if npix == 1*u.pixel:
            # Calculate new value based on area of FWHM and pixelscale
            # IRAF `ccdtime` uses this definition
    #        npix = 1.4*(fwhm / pixel_scale)**2
            npix = np.pi*(fwhm / pixel_scale)**2
            npix = int(round(npix.value)) * u.pixel

        print("npix=", npix)
        sobj2 = countrate * exp_time
        sky_countrate = sky_countrate * exp_time/obs_filter.equivwidth()
        dark_count = darkcurrent_rate * exp_time
        noise_ccd = np.sqrt(npix * (dark_count + readnoise_sq))
        print("Noise(Star)=", np.sqrt(sobj2))
        print("Noise( CCD)=", noise_ccd)
        ssky2 = background_rate * exp_time
        peak = sobj2/1.14/fwhm/fwhm*pixel_scale*pixel_scale


        print('Detected photons            from object =', sobj2, \
             ' (max',peak,' /pixel)')
        print('Detected photons/A/arcsec^2 from sky    =', sky_countrate)
        print('Detected photons/pixel      from sky    =', ssky2)
        snr = self._compute_snr(sobj2, ssky2, npix, dark_count, readnoise_sq)
        print('Signal-to-noise                         =', snr.value)
        return snr.value

    def _compute_snr(self, sobj2, ssky2, npix, dark_count, readnoise_sq):
        """Compute the signal-to-noise via the CCD equation
        Ref: Howell (2006), page 76
        """
        snr = sobj2 / np.sqrt(sobj2+npix*(ssky2+dark_count+readnoise_sq))

        return snr

    def exptime_from_ccd_snr(self, snr, V_mag, filtername, npix=1 * u.pixel,
                         n_background=np.inf * u.pixel,
                         background_rate=0 * (u.ct / u.pixel / u.s),
                         sky_mag=None,
                         darkcurrent_rate=0 * (u.ct / u.pixel / u.s)):

        countrate = self.photons_from_source(V_mag, filtername)

        # make sure the gain is in the correct units
        gain = self.instrument.ccd_gain
        if gain.unit in (u.electron / u.adu, u.photon / u.adu):
            gain = gain.value * (u.ct / u.adu)
        elif gain.unit != u.ct / u.adu:
            raise u.UnitsError('gain must have units of (either '
                               'astropy.units.ct, astropy.units.electron, or '
                               'astropy.units.photon) / astropy.units.adu')

        # get the countrate from the synphot observation object
        countrate = countrate / gain
        print(countrate)
        # define counts to be in ADU, which are not technically convertible in
        # astropy.units:
        countrate = countrate.value * (u.ct / u.s)

        # necessary for units to work in countrate calculation:
        if not hasattr(snr, 'unit'):
            snr = snr * np.sqrt(1 * u.ct)
        readnoise = self._get_shotnoise(self.instrument.ccd_readnoise)
        gain_err = self._get_shotnoise(self.instrument.ccd_gain * self.instrument.adc_error)
        print(readnoise, gain_err)
        if sky_mag is not None:
            # Compute countrate from given sky magnitude.
            # XXX need to refactor photons_from_source() to do this also
            sky_filtername = filtername
            if '::' in filtername:
                sky_filtername = filtername.split('::')[1]
            sky = self.site.sky_spectrum(sky_filtername)
            print("  Original sky=" ,sky(sky.avgwave()))
            sky2 = sky.normalize(sky_mag*units.VEGAMAG, self._V_band, vegaspec=self._vega)
            print("Normalized sky=", sky2(sky2.avgwave()), sky(sky.avgwave())*10**(-sky_mag/2.5))
            self._create_combined()
            waves, thru = self.combined._get_arrays(None)
            obs_filter = self.instrument.filterset[filtername]
            filter_waves, filter_trans = obs_filter._get_arrays(waves)
            signal_combined = self.combined/self.instrument.ccd_qe
            signal_waves, signal_thru = signal_combined._get_arrays(None)
            print('thruput (Atm*tel*instr)=', signal_thru.max())
            print('thruput (Atm*tel*instr*CCD)=', thru.max(), self.telescope.area.to(units.AREA)*thru.max())
            spec_elements = SpectralElement(Empirical1D, points=filter_waves, lookup_table = thru * filter_trans)

            # get the synphot observation object
            sky_obs = Observation(sky2, spec_elements, force='taper')

            sky_countrate = sky_obs.countrate(area=self.telescope.area)
            # Actually a count rate per arcsec^2 as the original magnitude is
            # also per sq. arcsec
            sky_countrate /= u.arcsec**2
            print("Sky (photons/AA/arcsec^2)=", sky_countrate, sky_countrate/(1.0*860)*u.AA)
            sky_per_pixel = sky_countrate * self.instrument.ccd_pixscale * self.instrument.ccd_pixscale
            sky_per_pixel /= u.pixel
            print("Sky (photons/pixel)=", sky_per_pixel)
            background_rate = sky_per_pixel

        npix = np.pi*(self.instrument.fwhm / self.instrument.ccd_pixscale)**2
        npix *= u.pixel
        print("npix=", npix)

        # solve t with the quadratic equation (pg. 77 of Howell 2006)
        A = countrate ** 2
        B = (-1) * snr ** 2 * (countrate + npix * (background_rate +
                                                   darkcurrent_rate))
        C = (-1) * snr ** 2 * npix * readnoise ** 2
        t = (-B + np.sqrt(B ** 2 - 4 * A * C)) / (2 * A)

        if gain_err.value > 1 or np.isfinite(n_background.value):
            from scipy.optimize import fsolve
            # solve t numerically
            t = fsolve(_t_with_small_errs, t, args=(background_rate,
                                                    darkcurrent_rate,
                                                    gain_err, readnoise, countrate,
                                                    npix, n_background))
            t = float(t) * u.s

        return t


    def _get_shotnoise(self, detector_property):
        """
        Returns the shot noise (i.e. non-Poissonion noise) in the correct
        units. ``detector_property`` must be a Quantity.
        """
        # Ensure detector_property is in the correct units:
        if detector_property.unit in (u.electron / u.pix,):
            detector_property = (detector_property.value/self.instrument.ccd_gain.value) * (u.ct/u.pix)
        detector_property = detector_property.to(u.ct / u.pixel)
        return detector_property.value * np.sqrt(1 * (u.ct / u.pixel))


    def _t_with_small_errs(self, t, background_rate, darkcurrent_rate, gain_err,
                           readnoise, countrate, npix, n_background):
        """
        Returns the full expression for the exposure time including the
        contribution to the noise from the background and the gain.
        """
        if not hasattr(t, 'unit'):
            t = t * u.s

        detector_noise = (background_rate * t + darkcurrent_rate * t +
                          gain_err ** 2 + readnoise ** 2)
        radicand = countrate * t + (npix * (1 + npix / n_background) *
                                    detector_noise)

        return countrate * t / np.sqrt(radicand)

    def _do_plot(self, waves, thru, filterlist=[], filterset=None, title='', left=300*u.nm, right=1200*u.nm, bottom=0.0, top=None):
        """Plot worker.

        Parameters
        ----------
        waves, thru : `~astropy.units.quantity.Quantity`
            Wavelength and throughput to plot.

        filterset : list
            List of filter bandpasses to plot

        kwargs :

        title : str
            Plot title.

        left, right : `None`, number or `Quantity`
            Minimum and maximum wavelengths to plot.
            If `None`, uses the whole range. If a number is given,
            must be in Angstrom. If an `~astropy.units.quantity.Quantity`
            is given, it will be converted tp the internal wave units.

        bottom, top : `None` or number
            Minimum and maximum flux/throughput to plot.
            If `None`, uses the whole range. If a number is given,
            it assumed to be in throughput (0..1).
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print('No matplotlib installation found; plotting disabled as a result.')
            return

        fig, ax = plt.subplots()
        ax.plot(waves, thru)

        # Plot filters
        for filtername in filterlist:
            filter_wave, filter_trans = filterset[filtername]._get_arrays(waves)
            filt_plot = ax.plot(filter_wave.to(self._internal_wave_unit), filter_trans*thru, linestyle='--', linewidth=1, label=filtername)

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
        ax.legend()

        plt.draw()

    def plot(self, filterlist=[], **kwargs):
        """Plot combined system

        Parameters
        ----------
        waves, thru : `~astropy.units.quantity.Quantity`
            Wavelength and throughput to plot.

        filterlist : list or str
            Filter name(s) to plot

        kwargs
            See :func:`do_plot`.

        """
        self._create_combined()
        waves, trans = self.combined._get_arrays(None)
        filterset = None
        # Handle single filter string case
        if isinstance(filterlist, str):
            filterlist = [filterlist,]
        if len(filterlist) > 0 and set(filterlist).issubset(set(self.instrument.filterlist)):
            filterset = self.instrument.filterset
        self._do_plot(waves.to(self._internal_wave_unit), trans, filterlist, filterset, **kwargs)

    def _create_combined(self):
        if self.combined is None:
            self.combined = self.site.transmission * self.telescope.reflectivity * self.instrument.transmission * self.instrument.ccd_qe

    def __repr__(self):
        prefixstr = '<' + self.__class__.__name__ + ' '
        compstr = "->".join([repr(x) for x in self.components])
        return f'{prefixstr}{compstr}>'
