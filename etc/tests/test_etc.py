import pytest
import os
import warnings
from copy import deepcopy

from astropy import units as u
warnings.simplefilter("ignore", pytest.PytestUnknownMarkWarning)
from astropy.tests.helper import assert_quantity_allclose
from synphot import units
from synphot.spectrum import SpectralElement, BaseUnitlessSpectrum

from etc import etc
from etc.models import ETCError

class TestETC:
    def test_initialize(self):
        test_etc = etc.ETC()

        assert '<ETC [Site(Maui, HI)]->[Telescope(FTN)]->[Instrument(FLOYDS)]>' == repr(test_etc)

    def test_initialize_specified_components(self):
        test_etc = etc.ETC(components=['sky', 'earth', 'underground cavern'])

        assert "<ETC 'sky'->'earth'->'underground cavern'>" == repr(test_etc)

    def test_initialize_from_dict(self):

        config_dict = { 'site' : {  'name': 'La Palma, Spain',
                                    'transmission' : 0.8709635899560807
                                 },
                        'telescope' : { 'name' : 'WHT 4.2m',
                                        'size': 4.2,
                                        'area': 12.47,
                                        'num_mirrors': 1,
                                        'reflectivity': 0.85},
                        'instrument' : { 'name' : 'WHT PFIP',
                                         'num_ar_coatings': 0,
                                         'num_inst_mirrors': 0,
                                         'num_inst_lenses': 1,
                                         'inst_lens_trans': 0.7,
                                         'fwhm': 1.0,
                                         'focal_scale': 17.5,
                                         'focal_scale_units': 'arcsec/mm',
                                         'ccd_qe': 0.8,
                                         'ccd_readnoise': 4.0,
                                         'ccd_gain': 0.9,
                                         'ccd_pixsize': 13.5
                                        }
                        }

        test_etc = etc.ETC(config_dict)

        assert "<ETC [Site(La Palma, Spain)]->[Telescope(WHT 4.2m)]->[Instrument(WHT PFIP)]>" == repr(test_etc)

    def test_initialize_config_file(self):
         test_config_file = os.path.abspath(os.path.join(__package__, 'etc', "tests", "data", "test1.toml"))
         test_etc = etc.ETC(config_file=test_config_file)
         assert '<ETC [Site(BPL)]->[Telescope(BPL 1-m)]->[Instrument(Sinistro)]>' == repr(test_etc)

         assert test_etc.site.name == "BPL"
         assert test_etc.telescope.name == "BPL 1-m"
         assert test_etc.instrument.name == "Sinistro"

    def test_properties(self):
        test_etc = etc.ETC()

        assert test_etc.site.name == "Maui, HI"
        assert test_etc.site.altitude == 3065 * u.m
        assert test_etc.telescope.name == "FTN"
        assert test_etc.instrument.name == "FLOYDS"

    def test_throughput_for_filter_noatmos_singlechannel(self):
        #                     Mirrors   AR       Lens
        expected_throughput = 0.85*0.85*0.995**2*0.93**3*0.9
        config_dict = {
                        'site' : {  'name': 'Tenerife, Spain',
                                    'transmission' : 0.9
                                 },
                        'telescope' : { 'name' : 'LCO 1.0m',
                                        'size': 1.0,
                                        'area': 0.63,
                                        'num_mirrors': 2,
                                        'reflectivity': 0.85},
                        'instrument' : {'name': 'Sinistro',
                                                       'inst_type': 'Imager',
                                                       'num_inst_mirrors': 0,
                                                       'num_inst_lenses': 3,
                                                       'inst_ar_coating_refl': 0.995,
                                                       'fwhm': 1.17,
                                                       'fwhm_units': 'arcsec',
                                                       'focal_scale': 25.9,
                                                       'focal_scale_units': 'arcsec/mm',
                                                       'filterlist': ['U', 'B', 'V', 'R', 'I', 'up', 'gp', 'rp', 'ip', 'zs'],
                                                       'ccd_qe': 0.9,
                                                       'ccd_xpixels': 4096,
                                                       'ccd_ypixels': 4096,
                                                       'ccd_pixsize': 15.0,
                                                       }
                                        }

        test_etc = etc.ETC(config_dict)

        assert isinstance(test_etc._throughput_for_filter('rp', atmos=False), BaseUnitlessSpectrum)
        throughput = test_etc._throughput_for_filter('rp', atmos=False)(550*u.nm)
        assert_quantity_allclose(expected_throughput, throughput)

    def test_throughput_for_filter_singlechannel(self):
        #                     Atmos Mirrors   AR       Lens
        expected_throughput = 0.9 * 0.85*0.85*0.995**2*0.93**3*0.9
        config_dict = {
                        'site' : {  'name': 'Tenerife, Spain',
                                    'transmission' : 0.9
                                 },
                        'telescope' : { 'name' : 'LCO 1.0m',
                                        'size': 1.0,
                                        'area': 0.63,
                                        'num_mirrors': 2,
                                        'reflectivity': 0.85},
                        'instrument' : {'name': 'Sinistro',
                                                       'inst_type': 'Imager',
                                                       'num_inst_mirrors': 0,
                                                       'num_inst_lenses': 3,
                                                       'inst_ar_coating_refl': 0.995,
                                                       'fwhm': 1.17,
                                                       'fwhm_units': 'arcsec',
                                                       'focal_scale': 25.9,
                                                       'focal_scale_units': 'arcsec/mm',
                                                       'filterlist': ['U', 'B', 'V', 'R', 'I', 'up', 'gp', 'rp', 'ip', 'zs'],
                                                       'ccd_qe': 0.9,
                                                       'ccd_xpixels': 4096,
                                                       'ccd_ypixels': 4096,
                                                       'ccd_pixsize': 15.0,
                                                       }
                                        }

        test_etc = etc.ETC(config_dict)

        assert isinstance(test_etc._throughput_for_filter('rp'), BaseUnitlessSpectrum)
        throughput = test_etc._throughput_for_filter('rp')(550*u.nm)
        assert_quantity_allclose(expected_throughput, throughput)

    def test_throughput_for_filter_dualchannel(self):
        #                        Atmos Mirrors   AR       Lens    CCD
        expected_gp_throughput = 0.9 * 0.85*0.85*0.995**2*0.93**3*0.85
        expected_ip_throughput = 0.9 * 0.85*0.85*0.995**2*0.93**3*0.8
        config_dict = {
                        'site' : {  'name': 'Tenerife, Spain',
                                    'transmission' : 0.9
                                 },
                        'telescope' : { 'name' : 'LCO 1.0m',
                                        'size': 1.0,
                                        'area': 0.63,
                                        'num_mirrors': 2,
                                        'reflectivity': 0.85},
                        'instrument' : {'name': 'DBI',
                                        'inst_type': 'Imager',
                                        'num_inst_mirrors': 0,
                                        'num_ar_coatings' : 2,
                                        'num_inst_lenses': 3,
                                        'inst_ar_coating_refl': 0.995,
                                        'fwhm': 1.17,
                                        'fwhm_units': 'arcsec',
                                        'focal_scale': 25.9,
                                        'focal_scale_units': 'arcsec/mm',
                                        'ccd_xpixels': 4096,
                                        'ccd_ypixels': 4096,
                                        'ccd_pixsize': 15.0,
                                        'channels' : { 'blue_channel' : {'filterlist': ['gp',], 'ccd_qe': 0.85},
                                                       'red_channel'  : {'filterlist': ['ip',], 'ccd_qe': 0.8}
                                                     }
                                       }
                        }

        test_etc = etc.ETC(config_dict)

        assert isinstance(test_etc._throughput_for_filter('gp'), BaseUnitlessSpectrum)
        throughput = test_etc._throughput_for_filter('gp')(500*u.nm)
        assert_quantity_allclose(expected_gp_throughput, throughput)

        assert isinstance(test_etc._throughput_for_filter('ip'), BaseUnitlessSpectrum)
        throughput = test_etc._throughput_for_filter('ip')(750*u.nm)
        assert_quantity_allclose(expected_ip_throughput, throughput)

    def test_channel_for_filter_singlechannel(self):
        config_dict = {'instrument' : {'name': 'Sinistro',
                                       'inst_type': 'Imager',
                                       'num_inst_mirrors': 0,
                                       'num_inst_lenses': 3,
                                       'inst_ar_coating_refl': 0.995,
                                       'fwhm': 1.17,
                                       'fwhm_units': 'arcsec',
                                       'focal_scale': 25.9,
                                       'focal_scale_units': 'arcsec/mm',
                                       'filterlist': ['U', 'B', 'V', 'R', 'I', 'up', 'gp', 'rp', 'ip', 'zs'],
                                       'ccd_qe': 0.9,
                                       'ccd_xpixels': 4096,
                                       'ccd_ypixels': 4096,
                                       'ccd_pixsize': 15.0,
                                       }
                        }

        test_etc = etc.ETC(config_dict)

        assert test_etc.instrument.name == "Sinistro"
        assert test_etc.instrument.filterlist == config_dict['instrument']['filterlist']

        assert test_etc._channel_for_filter('V') == 0

    def test_channel_for_filter_dualchannel(self):
        config_dict = {'instrument' : {'name': 'DBI',
                                       'inst_type': 'Imager',
                                       'num_inst_mirrors': 0,
                                       'num_inst_lenses': 3,
                                       'inst_ar_coating_refl': 0.995,
                                       'fwhm': 1.17,
                                       'fwhm_units': 'arcsec',
                                       'focal_scale': 25.9,
                                       'focal_scale_units': 'arcsec/mm',
                                       'ccd_xpixels': 4096,
                                       'ccd_ypixels': 4096,
                                       'ccd_pixsize': 15.0,
                                       'channels' : { 'blue_channel' : {'filterlist': ['gp'], 'ccd_qe': 0.85},
                                                      'red_channel'  : {'filterlist': ['ip'], 'ccd_qe': 0.8}
                                                    }
                                       }
                        }

        test_etc = etc.ETC(config_dict)

        assert test_etc.instrument.name == "DBI"
        assert test_etc.instrument.filterlist == ['gp', 'ip']

        assert test_etc._channel_for_filter('gp') == 0
        assert test_etc._channel_for_filter('ip') == 1


class TestComputeSNR:
    @pytest.fixture(autouse=True)
    def setUp(self):
        self.gain = 5.0*(u.electron/u.adu)
        self.readnoise = 5
        self.darkcurrent_rate = 22*(u.electron/u.pixel/u.hour)
        self.exp_time = 300*u.s

    def test_howell_example(self):
        expected_snr = 342.07825332

        test_etc = etc.ETC()

        signal_obj = 24013 * self.gain.value
        signal_sky = 620 * self.gain.value
        npix = 1
        signal_dark = self.darkcurrent_rate * self.exp_time
        signal_dark = signal_dark.decompose().value
        readnoise_sq = self.readnoise * self.readnoise

        snr = test_etc._compute_snr(signal_obj, signal_sky, npix, signal_dark, readnoise_sq)

        assert_quantity_allclose(expected_snr, snr, rtol=1e-5)


class TestPhotonsFromSource:
    @pytest.fixture(autouse=True)
    def setUp(self):

        self.test_config_file = os.path.abspath(os.path.join(__package__, 'etc', "tests", "data", "test1.toml"))
        self.WHT_config = { 'site' : {  'name': 'La Palma, Spain',
                                        'transmission' : 0.8709635899560807
                                     },
                            'telescope' : { 'name' : 'WHT 4.2m',
                                            'size': 4.2,
                                            'area': 12.47,
                                            'num_mirrors': 1,
                                            'reflectivity': 0.85},
                            'instrument' : { 'name' : 'WHT PFIP',
                                             'num_ar_coatings': 0,
                                             'num_inst_mirrors': 0,
                                             'num_inst_lenses': 1,
                                             'inst_lens_trans': 0.7,
                                             'filterlist' : ['WHT::V',],
                                             'fwhm': 1.0,
                                             'focal_scale': 17.5,
                                             'focal_scale_units': 'arcsec/mm',
                                             'ccd_qe': 0.8,
                                             'ccd_readnoise': 4.0,
                                             'ccd_gain': 0.9,
                                             'ccd_pixsize': 13.5
                                            }
                            }

        self.flux_mag0_V = 998.8095092773 * (u.ct/u.s/u.cm**2/u.AA)

    def test_WHT_V20(self):
        # For V=20
        V_mag = 20.0

        # Countrate/s from source is calculated as follows:
        # 10.**(-mag/2.5)
        #   * photons/sec/A/cm**2 giving mag = 0 at top of atm.
        #   * Atm transmission [10.**(-extin/2.5*airmass) ]
        #   * measured throughput of telescope/instrument
        #   * Unobstructed area of mirror (converted to cm**2)
        #   * effective bandwidth
        #   * quantum efficiency of detector
        expected_countrate = 10**(-V_mag/2.5) \
                             * self.flux_mag0_V \
                             * self.WHT_config['site']['transmission'] \
                             * self.WHT_config['telescope']['reflectivity'] \
                             * self.WHT_config['instrument']['inst_lens_trans'] \
                             * (self.WHT_config['telescope']['area']*u.m**2).to(u.cm**2) \
                             * 855.179810*u.AA \
                             * self.WHT_config['instrument']['ccd_qe'] \

        test_etc = etc.ETC(self.WHT_config)

        countrate = test_etc.photons_from_source(V_mag, 'V', 'WHT::V')

        # Can't get a very exact match as we are using as actual filter transmission
        # not exactly represented by the equivalent width above but it's within 0.27%
        assert_quantity_allclose(expected_countrate, countrate, rtol=3e-3)

    def test_bad_sourcespec(self):

        test_etc = etc.ETC(self.WHT_config)

        with pytest.raises(Exception) as execinfo:
            countrate = test_etc.photons_from_source(15, 'V', 'WHT::V', 3631*u.Jy)

        assert execinfo.type == ETCError
        assert execinfo.value.args[0] == "Invalid sourcespec; must be a SourceSpectrum"

    def test_sky_spectrum_V(self):

        # Value from SIGNAL using default V band m_0=3640 Jy but effective bandwidth=855.18 AA
        # from import of SVO WHT/PFIP.Har_V filter profile as a SpectralElement and equivwidth()
        expected_countrate = 76.738514 * (u.ct/u.s)

        test_etc = etc.ETC(self.WHT_config)

        sky_spec = test_etc.site.sky_spectrum('V')
        countrate = test_etc.photons_from_source(21.9, 'V', 'WHT::V', sky_spec)

        assert_quantity_allclose(expected_countrate, countrate, rtol=9e-04)

    def test_sky_spectrum_B(self):

        WHT_config = deepcopy(self.WHT_config)
        WHT_config['site']['transmission'] = 0.7943282347242815
        WHT_config['instrument']['filterlist'] = ['WHT::B']
        WHT_config['instrument']['ccd_qe'] = 0.82
        WHT_config['instrument']['inst_lens_trans'] = 0.8 * 0.9

        print(WHT_config)
        test_etc = etc.ETC(WHT_config)

        # From SIGNAL with W=704.17948 (from integration of SVO WHT::B filter)
        # and more accurate Planck constant (6.62607015d0 vs 6.6d0)
        expected_countrate = 41.5181084 * (u.ct/u.s)

        sky_spec = test_etc.site.sky_spectrum('B')
        countrate = test_etc.photons_from_source(22.7, 'B', 'WHT::B', sky_spec)

        # Difference is ~1.55% from differences in normalization with exact
        # Vega spectrum vs standard 10**(-0.4*skymag)
        assert_quantity_allclose(expected_countrate, countrate, rtol=2e-2)


class TestCCDSNR:
    @pytest.fixture(autouse=True)
    def setUp(self):

        self.test_config_file = os.path.abspath(os.path.join(__package__, 'etc', "tests", "data", "test1.toml"))
        self.WHT_config = { 'site' : {  'name': 'La Palma, Spain',
                                        'transmission' : 0.8709635899560807
                                     },
                            'telescope' : { 'name' : 'WHT 4.2m',
                                            'size': 4.2,
                                            'area': 12.47,
                                            'num_mirrors': 1,
                                            'reflectivity': 0.85},
                            'instrument' : { 'name' : 'WHT PFIP',
                                             'num_ar_coatings': 0,
                                             'num_inst_mirrors': 0,
                                             'num_inst_lenses': 1,
                                             'inst_lens_trans': 0.7,
                                             'filterlist' : ['WHT::V',],
                                             'fwhm': 1.0,
                                             'focal_scale': 17.5,
                                             'focal_scale_units': 'arcsec/mm',
                                             'ccd_qe': 0.8,
                                             'ccd_readnoise': 4.0,
                                             'ccd_gain': 0.9,
                                             'ccd_pixsize': 13.5
                                            }
                            }

    def test_LCO_1m_1s_V15_default_sky(self):
        expected_snr = 10.4210644

        test_etc = etc.ETC(self.test_config_file)

        snr = test_etc.ccd_snr(1, 15, 'V')

        assert_quantity_allclose(expected_snr, snr, rtol=1e-5)

    def test_LCO_1m_1s_bad_filter(self):

        test_etc = etc.ETC(self.test_config_file)

        with pytest.raises(Exception) as execinfo:
            snr = test_etc.ccd_snr(1, 15, 'Vibble')

        assert execinfo.type == ETCError
        assert execinfo.value.args[0] == "Filter name Vibble is invalid."

    def test_LCO_1m_1s_V15(self):
        expected_snr = 10.4210644

        test_etc = etc.ETC(self.test_config_file)

        snr = test_etc.ccd_snr(1, 15, 'V', sky_mag=21.8)

        assert_quantity_allclose(expected_snr, snr, rtol=1e-5)

    def test_LCO_1m_1s_R15(self):
        expected_snr = 10.815680017464334

        test_etc = etc.ETC(self.test_config_file)

        snr = test_etc.ccd_snr(1, 15, 'R', sky_mag=21.0)

        assert_quantity_allclose(expected_snr, snr, rtol=1e-5)

    def test_WHT_1s_B15(self):

        # Note Harris B filter used in SIGNAL is very much different from a
        # regular Johnson/Bessell B filter and the flux (in Jy) for mag=0
        # is also very different from that in SVO for the same filter
        # (4086 Jy vs 3856.92) see:
        # http://svo2.cab.inta-csic.es/theory/fps/index.php?id=WHT/PFIP.Har_B&&mode=browse&gname=WHT&gname2=PFIP#filter
        # Also very different to the 3925.77 Jy given by SVO for a Generic Bessel B see:
        # http://svo2.cab.inta-csic.es/theory/fps/index.php?id=Generic/Bessell.B&&mode=browse&gname=Generic&gname2=Bessell#filter
        # Value below from a recompile using 3857.0 Jy
        expected_snr = 203.09965515

        extin = 0.25 # mag/airmass (from SIGNAL)
        airmass = 1.5
        self.WHT_config['site']['transmission'] = 10.**(-0.4*extin*airmass)
        self.WHT_config['instrument']['filterlist'] = ['WHT::B',]
        self.WHT_config['instrument']['inst_lens_trans'] = 0.80 * 0.9
        self.WHT_config['instrument']['ccd_qe'] = 0.82

        test_etc = etc.ETC(self.WHT_config)

        snr = test_etc.ccd_snr(1, 15, 'WHT::B', sky_mag=22.7)

        assert_quantity_allclose(self.WHT_config['site']['transmission'], 0.707945764)
        assert_quantity_allclose(expected_snr, snr, rtol=3e-3)

    def test_WHT_5s_R15(self):
        # Also a mismatch for R although not as severe. Default is 3080 Jy;
        # Value below from a recompile using 3857.0 Jy
        expected_snr = 476.35858154

        extin = 0.09 # mag/airmass (from SIGNAL)
        airmass = 1.2
        self.WHT_config['site']['transmission'] = 10.**(-0.4*extin*airmass)
        self.WHT_config['instrument']['filterlist'] = ['WHT::R',]
        self.WHT_config['instrument']['inst_lens_trans'] = 1.0 * 0.7
        self.WHT_config['instrument']['ccd_qe'] = 0.74

        test_etc = etc.ETC(self.WHT_config)

        snr = test_etc.ccd_snr(5, 15, 'WHT::R', sky_mag=21.0)

        assert_quantity_allclose(expected_snr, snr, rtol=3e-3)

    def test_WHT_V_specific_flux_rates(self):
        expected_snr = 168.05989075

        # For V=20, values from customized SIGNAL
        V_flux = 441.58419800 * (u.photon/u.s)
        sky_rate = 4.28308773 * (u.photon/u.pixel/u.s)

        test_etc = etc.ETC(self.WHT_config)

        snr = test_etc.ccd_snr(100, V_flux, 'WHT::V', background_rate=sky_rate)

    def test_WHT_R_specific_flux_rates(self):
        expected_snr = 168.05989075

        # For R=20, values from customized SIGNAL
        R_flux = 486.03009033 * (u.photon/u.s)
        sky_rate = 10.79959011 * (u.photon/u.pixel/u.s)
        extin = 0.09 # mag/airmass (from SIGNAL)
        airmass = 1.2

        self.WHT_config['site']['transmission'] = 10.**(-0.4*extin*airmass)
        self.WHT_config['instrument']['filterlist'] = ['WHT::R',]
        self.WHT_config['instrument']['ccd_qe'] = 0.74

        test_etc = etc.ETC(self.WHT_config)

        snr = test_etc.ccd_snr(100, R_flux, 'WHT::R', background_rate=sky_rate)
