import pytest
import os
import warnings

from astropy import units as u
warnings.simplefilter("ignore", pytest.PytestUnknownMarkWarning)
from astropy.tests.helper import assert_quantity_allclose

from etc import etc


def test_initialize():
    test_etc = etc.ETC()

    assert '<ETC [Site(Maui, HI)]->[Telescope(FTN)]->[Instrument(FLOYDS)]>' == repr(test_etc)

def test_initialize_specified_components():
    test_etc = etc.ETC(components=['sky', 'earth', 'underground cavern'])

    assert "<ETC 'sky'->'earth'->'underground cavern'>" == repr(test_etc)

def test_initialize_config_file():
     test_config_file = os.path.abspath(os.path.join(__package__, 'etc', "tests", "data", "test1.toml"))
     test_etc = etc.ETC(config_file=test_config_file)
     assert '<ETC [Site(BPL)]->[Telescope(BPL 1-m)]->[Instrument(Sinistro)]>' == repr(test_etc)

     assert test_etc.site.name == "BPL"
     assert test_etc.telescope.name == "BPL 1-m"
     assert test_etc.instrument.name == "Sinistro"

def test_properties():
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

class TestCCDSNR:
    @pytest.fixture(autouse=True)
    def setUp(self):

        self.test_config_file = os.path.abspath(os.path.join(__package__, 'etc', "tests", "data", "test1.toml"))

    def test_LCO_1m_1s_V15(self):
        expected_snr = 10.4210644

        test_etc = etc.ETC(self.test_config_file)

        snr = test_etc.ccd_snr(1, 15, 'V', sky_mag=21.8)

        assert_quantity_allclose(expected_snr, snr, rtol=1e-5)
