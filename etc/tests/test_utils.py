import pytest
import os
import warnings

import toml
from astropy import units as u
warnings.simplefilter("ignore", pytest.PytestUnknownMarkWarning)
from astropy.tests.helper import assert_quantity_allclose
from synphot.spectrum import SpectralElement, BaseUnitlessSpectrum, SourceSpectrum
from synphot import units

from etc.utils import read_element, read_eso_spectra

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

    def test_read_element_angstroms(self):
        test_fp = os.path.abspath(os.path.join(__package__, 'etc', "tests", "data", "test_atmos_angstroms.dat"))

        element = read_element(test_fp)

        assert_quantity_allclose(element.waveset[0], 3200 * u.AA)
        assert_quantity_allclose(element(element.waveset[0]), 0.306, rtol=1e-3)

    def test_read_element_as_sourcespectrum(self):
        test_fp = os.path.abspath(os.path.join(__package__, 'etc', "tests", "data", "test_radiance.dat"))

        element = read_element(test_fp, element_type='spectrum')

        assert isinstance(element, SourceSpectrum)
        assert element.waveset[0].unit == u.AA
        assert_quantity_allclose(element.waveset[0], 300 * u.nm)
        assert_quantity_allclose(element(element.waveset[0]), 16.4805*units.PHOTLAM, rtol=1e-3)

    def test_read_element_as_sourcespectrum_ownfluxunits(self):
        test_fp = os.path.abspath(os.path.join(__package__, 'etc', "tests", "data", "test_radiance.dat"))

        eso_unit = units.u.photon/units.u.s/units.u.m**2/units.u.um

        element = read_element(test_fp, element_type='spectrum', flux_units=eso_unit)

        assert isinstance(element, SourceSpectrum)
        assert element.waveset[0].unit == u.AA
        assert element(element.waveset[0]).unit == units.PHOTLAM
        assert_quantity_allclose(element.waveset[0], 300 * u.nm)
        assert_quantity_allclose(element(element.waveset[0]), 16.4805*eso_unit, rtol=1e-3)


class TestReadESOSpectra():

    @classmethod
    def setup_class(cls):

        cls.test_fits_fp = os.path.abspath(os.path.join(__package__, 'etc', "tests", "data", "test_spectrum.fits"))
        cls.funit = u.erg / (u.cm**2 * u.s * u.AA)

    def test_spectra(self):

        source_spec = read_eso_spectra(self.test_fits_fp)

        assert type(source_spec) == SourceSpectrum
        assert_quantity_allclose(source_spec.waveset[0], 3451.65 * u.AA)
        assert_quantity_allclose(source_spec(source_spec.waveset[0], flux_unit='FLAM'), 4.341547e-16*self.funit)
