import pytest
import os
import warnings

import toml
from astropy import units as u
warnings.simplefilter("ignore", pytest.PytestUnknownMarkWarning)
from astropy.tests.helper import assert_quantity_allclose
from synphot.spectrum import SpectralElement, BaseUnitlessSpectrum

from etc.utils import read_element

class TestReadElement:

    def test_read_element(self):
        test_fp = os.path.abspath(os.path.join(__package__, 'etc', "tests", "data", "test_mirror.dat"))

        element = read_element(test_fp)

        assert_quantity_allclose(element.waveset[0], 2750 * u.AA)
        assert_quantity_allclose(element(element.waveset[0]), 0.64)

    def test_read_element_microns(self):
        test_fp = os.path.abspath(os.path.join(__package__, 'etc', "tests", "data", "test_mirror_microns.dat"))

        element = read_element(test_fp)

        assert_quantity_allclose(element.waveset[0], 3000 * u.AA)
        assert_quantity_allclose(element(element.waveset[0]), 0.751)

    def test_read_element_microns_passed_waveunit(self):
        test_fp = os.path.abspath(os.path.join(__package__, 'etc', "tests", "data", "test_mirror_microns.dat"))

        element = read_element(test_fp, wave_units=u.micron)

        assert_quantity_allclose(element.waveset[0], 3000 * u.AA)
        assert_quantity_allclose(element(element.waveset[0]), 0.751)

    def test_read_element_microns_passed_fluxunit(self):
        test_fp = os.path.abspath(os.path.join(__package__, 'etc', "tests", "data", "test_mirror_microns.dat"))

        element = read_element(test_fp, flux_units=u.percent)

        assert_quantity_allclose(element.waveset[0], 3000 * u.AA)
        assert_quantity_allclose(element(element.waveset[0]), 0.00751) # File is not in percent
