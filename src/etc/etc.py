import os
import sys
import toml
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

import numpy as np
from astropy import units as u
from astropy.utils.exceptions import AstropyUserWarning

from synphot import SourceSpectrum, SpectralElement, units
from synphot.models import Empirical1D
from synphot.observation import Observation

from . import data
from .models import Site, Telescope, Instrument
from .config import conf
from .utils import sptype_to_pickles_standard, ETCError


class ETC(object):
    _internal_wave_unit = u.nm
    _sun_raw = SourceSpectrum.from_file(os.path.expandvars(conf.sun_file))
    _vega = SourceSpectrum.from_file(os.path.expandvars(conf.vega_file))
    _V_band = SpectralElement.from_filter("johnson_v")

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
        self.combined_noatmos = None
        self.obs = None

    def _create_components_from_config(self, config):
        class_mapping = {"site": Site, "telescope": Telescope, "instrument": Instrument}
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

    def _convert_filtername(self, filtername):
        """Converts system+filter name (e.g. 'WHT::B' to a bare filter name
        which is returned"""

        new_filtername = filtername
        if "::" in filtername:
            new_filtername = filtername.split("::")[-1]
        return new_filtername

    def _map_filter_to_standard(self, filtername):
        """Maps a passed <filtername> to a standard Johnson/Cousins/Bessell
        filter profile. Returns a SpectralElement, defaulting to V"""

        band_mapping = {
            "U": "johnson_u",
            "B": "johnson_b",
            "V": "johnson_v",
            "R": conf.bessell_R_file,
            "Rc": "cousins_r",
            "I": conf.bessell_I_file,
            "Ic": "cousins_i",
            "J": "bessel_j",
            "H": "bessel_h",
            "K": "bessel_k",
            "G": conf.gaiadr2_g_file,
        }
        band_name = band_mapping.get(self._convert_filtername(filtername), "johnson_v")
        if band_name.startswith("http://svo"):
            band = SpectralElement.from_file(band_name)
        else:
            band = SpectralElement.from_filter(band_name)

        return band

    def pickles_to_source_spec(self, sp_type):
        """Returns a SourceSpectrum from the passed <sp_type> if it is found
        in the Pickles atlas, otherwise None is returned.
        The path to the library (which can be remote e.g.
        https://ssb.stsci.edu/trds/grid/pickles/dat_uvi) is given in
        `config.Conf.pickles_library_path"""

        source_spec = None
        filename = sptype_to_pickles_standard(sp_type)
        if filename:
            filepath = os.path.expandvars(os.path.join(conf.pickles_library_path, filename))
            source_spec = SourceSpectrum.from_file(filepath)
        return source_spec

    def sso_to_source_spec(self, tax_type):
        """Returns a SourceSpectrum from the passed <tax_type> if it is found
        in the Bus-DeMeo taxonomy, otherwise None is returned.
        The spectrum is produced by multiplying the reflectance spectra by
        a Kurucz model for the Sun so will need to be normalized"""

        source_spec = None

        if tax_type.lower().startswith("sso::") is False:
            tax_type = "sso::" + tax_type.strip()
        config_item = conf.source_mapping.get(tax_type, None)
        if config_item is not None:
            filename = config_item()
            file_path = os.path.expandvars(filename)
            if not os.path.exists(file_path):
                file_path = str(pkg_resources.files("etc.data").joinpath(filename))
            # Flux unit for LSST throughputs (*almost* FLAM but nm not Angstroms)
            lsst_funit = u.erg / u.cm**2 / u.s / u.nm
            source_spec = SourceSpectrum.from_file(
                file_path, wave_unit=u.nm, flux_unit=lsst_funit, header_start=1
            )
            source_spec.meta["header"]["source"] = config_item.description
            source_spec.meta["header"]["filename"] = filename

        return source_spec

    def photons_from_source(
        self, mag, mag_filter, filtername, source_spec=None, force_overlap="taper", normalize=True
    ):
        """Computes the photons coming from the source given by [source_spec] (Vega
        is assumed if not given) of the given <mag> in <mag_filter> (e.g. for "V=20",
        mag=20, mag_filter='V') when observed through the instrument's filter
        given by <filtername>
        """
        print("photons_from_source:", filtername)
        standard_filter = self._convert_filtername(filtername)
        band = self._map_filter_to_standard(mag_filter)
        print(band.meta.get("expr", "Unknown"), type(source_spec))
        if source_spec is None:
            source_spec = self._vega
        if type(source_spec) != SourceSpectrum:
            raise ETCError("Invalid sourcespec; must be a SourceSpectrum")

        # XXX Todo: allow support for AB mags here
        # test line below for comparison with SIGNAL. Difference in mag when
        # using proper normalization is ~mag-0.0155 (for B=22.7)
        #        source_spec_norm = source_spec*10**(-0.4*mag)
        if normalize is True:
            source_spec_norm = source_spec.normalize(mag * units.VEGAMAG, band, vegaspec=self._vega)
        else:
            source_spec_norm = source_spec

        self._create_combined()
        filter_waves, filter_trans = self.instrument.filterset[filtername]._get_arrays(None)
        throughput = self._throughput_for_filter(filtername)
        waves, thru = throughput._get_arrays(filter_waves)
        spec_elements = SpectralElement(Empirical1D, points=filter_waves, lookup_table=filter_trans * thru)

        # get the synphot observation object
        synphot_obs = Observation(source_spec_norm, spec_elements, force=force_overlap)
        if self.obs is None:
            self.obs = synphot_obs
        countrate = synphot_obs.countrate(area=self.telescope.area)

        return countrate

    def ccd_snr(
        self,
        exp_time,
        V_mag,
        filtername,
        source_spec=None,
        npix=1 * u.pixel,
        n_background=np.inf * u.pixel,
        background_rate=0 * (u.ph / u.pixel / u.s),
        sky_mag=None,
        darkcurrent_rate=0 * (u.ph / u.pixel / u.s),
        normalize=True,
    ):
        """Calculate the SNR in a given exposure time <exp_time> for a given <V_mag>"""

        if filtername not in self.instrument.filterlist:
            raise ETCError("Filter name {0} is invalid.".format(filtername))

        try:
            exp_time = exp_time.to(u.s)
        except AttributeError:
            exp_time = exp_time * u.s

        # Get countrate per second for this magnitude in the filter
        # Check first if we weren't given a rate directly
        try:
            countrate = V_mag.to(u.photon / u.s)
        except (AttributeError, u.UnitConversionError):
            countrate = self.photons_from_source(
                V_mag, "V", filtername, source_spec=source_spec, normalize=normalize
            )

        # define countrate to be in photons/sec, which can be treated as the
        # same as (photo)electrons for CCDs but are not technically convertible in
        # astropy.units:
        countrate = countrate.value * (u.photon / u.s)

        print("Counts/s=", countrate)

        readnoise = self.instrument.ccd_readnoise.value * (u.photon / u.pixel)
        readnoise_sq = readnoise.value**2 * (u.photon / u.pixel)
        gain_err = self._get_shotnoise(self.instrument.ccd_gain * self.instrument.adc_error)
        print("Readnoise=", readnoise, readnoise_sq, gain_err)

        fwhm = self.instrument.fwhm
        pixel_scale = self.instrument.ccd_pixscale
        pixel_area = pixel_scale * pixel_scale

        obs_filter = self.instrument.filterset[filtername]

        if background_rate > 0:
            print("Specific rate given")
            sky_countrate = background_rate / pixel_area / obs_filter.equivwidth()
        else:
            # Obtain sky value based on flux at mag=0 for given filter
            sky_filtername = self._convert_filtername(filtername)
            sky = self.site.sky_spectrum(sky_filtername)
            print("  Original sky for =", sky(sky.avgwave()), sky.meta.get("header", {}).get("filename", ""))
            prefix = "Using"
            sky2 = None
            if sky_mag is None:
                # Try and obtain sky magnitude from Site. This should probably
                # move into sky_spectrum()
                sky_mag = self.site.sky_mags.get(sky_filtername, None)
                prefix = "Determined new"
                if sky_mag is None:
                    if self.site.radiance is None:
                        raise ETCError(
                            "Could not determine a valid sky magnitude for filter: {0}".format(sky_filtername)
                        )
                    else:
                        print("Attempting synthesized sky magnitude from radiance spectrum")
                        normalizing_sky_mag = self.site.sky_mags.get("V", None)
                        V_band = self._map_filter_to_standard("V")
                        print("Normalizing sky spectrum to {} in V band".format(normalizing_sky_mag))
                        radiance_norm = sky.normalize(
                            normalizing_sky_mag * units.VEGAMAG, V_band, vegaspec=self._vega
                        )
                        obs_rad_norm = Observation(radiance_norm, obs_filter, force="taper")
                        sky_mag = obs_rad_norm.effstim(units.VEGAMAG, vegaspec=self._vega)
                        sky_mag = sky_mag.value
                        sky2 = radiance_norm

            print("{} sky magnitude of {:.2f} for {} band".format(prefix, sky_mag, sky_filtername))
            if sky2 is None:
                # Compute countrate from given sky magnitude.
                # XXX need to refactor photons_from_source() to do this also
                sky_filter = self._map_filter_to_standard(sky_filtername)
                print("sky_filter", sky_filter.meta.get("expr", "Unknown"), type(sky))

                sky2 = sky.normalize(sky_mag * units.VEGAMAG, sky_filter, vegaspec=self._vega)
            print("Normalized sky=", sky2(sky2.avgwave()), sky2(sky2.avgwave()) * 10 ** (-sky_mag / 2.5))
            self._create_combined()
            throughput = self._throughput_for_filter(filtername, atmos=False)
            waves, thru = throughput._get_arrays(None)
            filter_waves, filter_trans = obs_filter._get_arrays(waves)
            print("thruput=", thru.max(), self.telescope.area.to(units.AREA) * thru.max())
            spec_elements = SpectralElement(
                Empirical1D, points=filter_waves, lookup_table=thru * filter_trans
            )

            # get the synphot observation object
            warnings.simplefilter("ignore", category=AstropyUserWarning)
            sky_obs = Observation(sky2, spec_elements, force="taper")

            sky_countrate = sky_obs.countrate(area=self.telescope.area)
            # Actually a count rate per arcsec^2 as the original magnitude is
            # also per sq. arcsec. Also set to photons/s to match signal from object
            sky_countrate = sky_countrate.value * (u.photon / u.s / u.arcsec**2)
            print("Sky (photons/AA/arcsec^2)=", sky_countrate, sky_countrate / obs_filter.equivwidth())
            sky_per_pixel = sky_countrate * pixel_area
            sky_per_pixel /= u.pixel
            print("Sky (photons/pixel)=", sky_per_pixel)
            background_rate = sky_per_pixel

        # For spectrographs, calculate fraction of light passing through slit
        vign = self.instrument.slit_vignette()

        if npix == 1 * u.pixel:
            # Calculate new value based on area of FWHM and pixelscale
            # IRAF `ccdtime` uses this definition. Pass in a suitable
            # value for npix into this method to emulate this.
            #        npix = 1.4*(fwhm / pixel_scale)**2
            npix = np.pi * (fwhm / pixel_scale) ** 2
            print("npix (fractional pixels)=", npix)
            npix = int(round(npix.value)) * u.pixel

        print("npix (whole pixels)=", npix)
        sobj2 = countrate * exp_time
        sky_countrate = sky_countrate * exp_time / obs_filter.equivwidth()
        dark_count = darkcurrent_rate.to(u.photon / u.pix / u.s) * exp_time
        print("Noise(Star)=", np.sqrt(sobj2))
        print("Noise(Dark)=", np.sqrt(npix * dark_count), dark_count)
        print("Noise( CCD)=", np.sqrt(npix * readnoise_sq))
        noise_ccd = np.sqrt(npix * (dark_count + readnoise_sq))
        ssky2 = background_rate * exp_time
        peak = sobj2 / 1.14 / fwhm / fwhm * pixel_area

        print("Detected photons            from object =", sobj2, " (max", peak, " /pixel)")
        print("Detected photons/A/arcsec^2 from sky    =", sky_countrate)
        print("Detected photons/pixel      from sky    =", ssky2)
        snr = self._compute_snr(sobj2, ssky2, npix, dark_count, readnoise_sq)
        print("Signal-to-noise                         =", snr.value)
        return snr.value

    def _compute_snr(self, sobj2, ssky2, npix, dark_count, readnoise_sq):
        """Compute the signal-to-noise via the CCD equation
        Ref: Howell (2006), page 76
        """
        snr = sobj2 / np.sqrt(sobj2 + npix * (ssky2 + dark_count + readnoise_sq))

        return snr

    def exptime_from_ccd_snr(
        self,
        snr,
        V_mag,
        filtername,
        npix=1 * u.pixel,
        n_background=np.inf * u.pixel,
        background_rate=0 * (u.ct / u.pixel / u.s),
        sky_mag=None,
        darkcurrent_rate=0 * (u.ct / u.pixel / u.s),
    ):
        countrate = self.photons_from_source(V_mag, filtername)

        # make sure the gain is in the correct units
        gain = self.instrument.ccd_gain
        if gain.unit in (u.electron / u.adu, u.photon / u.adu):
            gain = gain.value * (u.ct / u.adu)
        elif gain.unit != u.ct / u.adu:
            raise u.UnitsError(
                "gain must have units of (either "
                "astropy.units.ct, astropy.units.electron, or "
                "astropy.units.photon) / astropy.units.adu"
            )

        # get the countrate from the synphot observation object
        countrate = countrate / gain
        print(countrate)
        # define counts to be in ADU, which are not technically convertible in
        # astropy.units:
        countrate = countrate.value * (u.ct / u.s)

        # necessary for units to work in countrate calculation:
        if not hasattr(snr, "unit"):
            snr = snr * np.sqrt(1 * u.ct)
        readnoise = self._get_shotnoise(self.instrument.ccd_readnoise)
        gain_err = self._get_shotnoise(self.instrument.ccd_gain * self.instrument.adc_error)
        print(readnoise, gain_err)
        if sky_mag is not None:
            # Compute countrate from given sky magnitude.
            # XXX need to refactor photons_from_source() to do this also
            sky_filtername = filtername
            if "::" in filtername:
                sky_filtername = filtername.split("::")[-1]
            sky = self.site.sky_spectrum(sky_filtername)
            print("  Original sky=", sky(sky.avgwave()))
            sky2 = sky.normalize(sky_mag * units.VEGAMAG, self._V_band, vegaspec=self._vega)
            print("Normalized sky=", sky2(sky2.avgwave()), sky(sky.avgwave()) * 10 ** (-sky_mag / 2.5))
            self._create_combined()
            waves, thru = self.combined._get_arrays(None)
            obs_filter = self.instrument.filterset[filtername]
            filter_waves, filter_trans = obs_filter._get_arrays(waves)
            signal_combined = self.combined / self.instrument.ccd_qe
            signal_waves, signal_thru = signal_combined._get_arrays(None)
            print("thruput (Atm*tel*instr)=", signal_thru.max())
            print("thruput (Atm*tel*instr*CCD)=", thru.max(), self.telescope.area.to(units.AREA) * thru.max())
            spec_elements = SpectralElement(
                Empirical1D, points=filter_waves, lookup_table=thru * filter_trans
            )

            # get the synphot observation object
            with warnings.simplefilter("ignore", category=AstropyUserWarning):
                sky_obs = Observation(sky2, spec_elements)

            sky_countrate = sky_obs.countrate(area=self.telescope.area)
            # Actually a count rate per arcsec^2 as the original magnitude is
            # also per sq. arcsec
            sky_countrate /= u.arcsec**2
            print("Sky (photons/AA/arcsec^2)=", sky_countrate, sky_countrate / (1.0 * 860) * u.AA)
            sky_per_pixel = sky_countrate * self.instrument.ccd_pixscale * self.instrument.ccd_pixscale
            sky_per_pixel /= u.pixel
            print("Sky (photons/pixel)=", sky_per_pixel)
            background_rate = sky_per_pixel

        npix = np.pi * (self.instrument.fwhm / self.instrument.ccd_pixscale) ** 2
        npix *= u.pixel
        print("npix=", npix)

        # solve t with the quadratic equation (pg. 77 of Howell 2006)
        A = countrate**2
        B = (-1) * snr**2 * (countrate + npix * (background_rate + darkcurrent_rate))
        C = (-1) * snr**2 * npix * readnoise**2
        t = (-B + np.sqrt(B**2 - 4 * A * C)) / (2 * A)

        if gain_err.value > 1 or np.isfinite(n_background.value):
            from scipy.optimize import fsolve

            # solve t numerically
            t = fsolve(
                _t_with_small_errs,
                t,
                args=(background_rate, darkcurrent_rate, gain_err, readnoise, countrate, npix, n_background),
            )
            t = float(t) * u.s

        return t

    def _get_shotnoise(self, detector_property):
        """
        Returns the shot noise (i.e. non-Poissonion noise) in the correct
        units. ``detector_property`` must be a Quantity.
        """
        # Ensure detector_property is in the correct units:
        if detector_property.unit in (u.electron / u.pix,):
            detector_property = (detector_property.value / self.instrument.ccd_gain.value) * (u.ct / u.pix)
        detector_property = detector_property.to(u.ct / u.pixel)
        return detector_property.value * np.sqrt(1 * (u.ct / u.pixel))

    def _t_with_small_errs(
        self, t, background_rate, darkcurrent_rate, gain_err, readnoise, countrate, npix, n_background
    ):
        """
        Returns the full expression for the exposure time including the
        contribution to the noise from the background and the gain.
        """
        if not hasattr(t, "unit"):
            t = t * u.s

        detector_noise = background_rate * t + darkcurrent_rate * t + gain_err**2 + readnoise**2
        radicand = countrate * t + (npix * (1 + npix / n_background) * detector_noise)

        return countrate * t / np.sqrt(radicand)

    def _do_plot(
        self,
        waves,
        thru,
        filterlist=[],
        filterset=None,
        title="",
        plot_label="",
        left=300 * u.nm,
        right=1200 * u.nm,
        bottom=0.0,
        top=None,
    ):
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

        plot_label : str
            Line plot label for legend

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
            print("No matplotlib installation found; plotting disabled as a result.")
            return

        fig, ax = plt.subplots()
        ax.plot(waves, thru, label=plot_label)

        # Plot filters
        for filtername in filterlist:
            filter_wave, filter_trans = filterset[filtername]._get_arrays(waves)
            filt_plot = ax.plot(
                filter_wave.to(self._internal_wave_unit),
                filter_trans * thru,
                linestyle="--",
                linewidth=1,
                label=filtername,
            )

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
        ax.set_xlabel("Wavelength ({0})".format(wave_unit))

        yu = thru.unit
        if yu is u.dimensionless_unscaled:
            ax.set_ylabel("Unitless")
        elif yu is u.percent:
            ax.set_ylabel("Efficiency ({0})".format(yu))
        else:
            ax.set_ylabel("Flux ({0})".format(yu))

        if title:
            ax.set_title(title)

        # Turn on minorticks on both axes
        ax.minorticks_on()
        if len(filterlist) > 0 or plot_label != "":
            ax.legend()

        plt.draw()

    def plot(self, filterlist=[], **kwargs):
        """Plot combined system throughput and filters. By default (if
        `filterlist=[]`) then all instrument filters are plotted. To plot
        specific filters, set `filterlist` to a string of the filtername or
        a list of filternames. To disable filter plotting completely, set
        `filterlist=None`

        Parameters
        ----------
        filterlist : list, str or None
            Filter name(s) to plot

        kwargs :
            atmos : bool
                Whether to include the atmosphere or not
            Also see :func:`do_plot`.

        """
        self._create_combined()
        if kwargs.get("atmos", True) is True:
            waves, trans = self.combined[0]._get_arrays(None)
        else:
            waves, trans = self.combined_noatmos[0]._get_arrays(None)
        for channel_num in range(1, self.instrument.num_channels):
            if kwargs.get("atmos", True) is True:
                chan_trans = self.combined[channel_num](waves)
            else:
                chan_trans = self.combined_noatmos[channel_num](waves)
            trans *= chan_trans

        if "atmos" in kwargs:
            del kwargs["atmos"]
        filterset = None
        # Handle single filter string case
        if isinstance(filterlist, str):
            filterlist = [
                filterlist,
            ]
        elif filterlist == []:
            filterlist = self.instrument.filterlist
        elif filterlist is None:
            filterlist = []
        if len(filterlist) > 0 and set(filterlist).issubset(set(self.instrument.filterlist)):
            filterset = self.instrument.filterset
        self._do_plot(waves.to(self._internal_wave_unit), trans, filterlist, filterset, **kwargs)

    def _efficiency_for_filter(self, filtername, atmos=True):
        """Calculate combined system throughput for filtername.

        Parameters
        ----------
        filtername : str
            Filter name to calculate efficiency for

        kwargs :
            atmos : bool
                Whether to include the atmosphere or not (default is True)

        Returns
        -------
        waves : Quantity
            Wavelength array (with units of _internal_wave_unit)
        thru : Quantity
            Combined system efficiency (units of %)
        """

        waves = None
        thru = None
        # Handle single filter string case
        if filtername in self.instrument.filterlist:
            throughput = self._throughput_for_filter(filtername, atmos)

            efficiency = throughput * self.instrument.filterset[filtername]
            # Add slit/fiber vignetting (if applicable)
            vignetting = self.instrument.slit_vignette()
            #            print("Vignetting", vignetting, (1.0-vignetting)*100., self.instrument.fiber_diameter, self.instrument.fwhm)
            efficiency = efficiency * vignetting

            waves = efficiency.waveset
            thru = efficiency(waves) * 100.0 * u.percent
            waves = waves.to(self._internal_wave_unit)
        else:
            raise ETCError("Filter name {0} is invalid.".format(filtername))

        return waves, thru

    def plot_efficiency(self, filtername, **kwargs):
        """Plot combined system throughput for filtername.

        Parameters
        ----------
        filtername : str
            Filter name to plot

        kwargs :
            atmos : bool
                Whether to include the atmosphere or not (default is True)
            Also see :func:`do_plot` for additional options.

        """

        # Handle single filter string case
        if filtername in self.instrument.filterlist:
            waves, thru = self._efficiency_for_filter(filtername, kwargs.get("atmos", True))

            if kwargs.get("atmos", True) is True:
                atmos_presence = "incl."
            else:
                atmos_presence = "excl."
            if "atmos" in kwargs:
                del kwargs["atmos"]
            if "title" not in kwargs:
                kwargs["title"] = "System Efficiency {} Atmosphere".format(atmos_presence)

            self._do_plot(waves, thru, filterlist=[], filterset=[], **kwargs)
        else:
            raise ETCError("Filter name {0} is invalid.".format(filtername))

    def _channel_for_filter(self, filtername):
        """Maps the passed <filtername> to Instrument channel number which is
        returned"""

        channel_name = self.instrument.filter2channel_map.get(filtername, "default")
        try:
            index = list(self.instrument.channelset).index(channel_name)
        except ValueError:
            print("Filter {} not in instrument channel map".format(filtername))
        return index

    def _throughput_for_filter(self, filtername, atmos=True):
        self._create_combined()
        if atmos is True:
            throughput = self.combined
        else:
            throughput = self.combined_noatmos
        # Map filtername to instrument channel and correct throughput
        index = self._channel_for_filter(filtername)
        throughput = throughput[index]

        return throughput

    def _create_combined(self):
        if self.combined_noatmos is None:
            self.combined_noatmos = []
            for channel in self.instrument.channelset.values():
                channel_combined = (
                    self.telescope.reflectivity
                    * self.instrument.transmission
                    * channel.transmission
                    * channel.ccd_qe
                )
                self.combined_noatmos.append(channel_combined)
        if self.combined is None:
            self.combined = []
            for channel in self.combined_noatmos:
                channel_combined = channel * self.site.transmission
                self.combined.append(channel_combined)
        # Collapse lists back down to single instance if only one channel
        # if self.instrument.num_channels == 1:
        # self.combined_noatmos = self.combined_noatmos[0]
        # self.combined = self.combined[0]

    def __repr__(self):
        prefixstr = "<" + self.__class__.__name__ + " "
        compstr = "->".join([repr(x) for x in self.components])
        return f"{prefixstr}{compstr}>"
