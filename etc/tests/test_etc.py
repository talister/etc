import pytest
import os
import warnings

from astropy import units as u
warnings.simplefilter("ignore", pytest.PytestUnknownMarkWarning)
from astropy.tests.helper import assert_quantity_allclose
from synphot import units

from etc import etc

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

        countrate = test_etc.photons_from_source(V_mag, 'WHT::V')

        # Can't get a very exact match as we are using as actual filter transmission
        # not exactly represented by the equivalent width above but it's within 0.27%
        assert_quantity_allclose(expected_countrate, countrate, rtol=3e-3)

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

    def test_LCO_1m_1s_V15(self):
        expected_snr = 10.4210644

        test_etc = etc.ETC(self.test_config_file)

        snr = test_etc.ccd_snr(1, 15, 'V', sky_mag=21.8)

        assert_quantity_allclose(expected_snr, snr, rtol=1e-5)

    def test_specific_flux_rates(self):
        expected_snr = 168.05989075

        # For V=20, values from customized SIGNAL
        V_flux = 441.58419800 * (u.photon/u.s)
        sky_rate = 4.28308773 * (u.photon/u.pixel/u.s)

        test_etc = etc.ETC(self.WHT_config)

        snr = test_etc.ccd_snr(100, V_flux, 'WHT::V', background_rate=sky_rate)
