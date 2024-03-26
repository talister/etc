import pytest
import os
import sys
import warnings

try:
    if sys.version_info[0] == 3 and sys.version_info[1] < 10:
        # Python 3.9 version of importlib.resources.files behaves differently
        from importlib_resources import files
    else:
        from importlib.resources import files
except ModuleNotFoundError:
    from importlib_resources import files

import toml
from astropy import units as u

warnings.simplefilter("ignore", pytest.PytestUnknownMarkWarning)
from astropy.tests.helper import assert_quantity_allclose
from synphot.spectrum import SpectralElement, BaseUnitlessSpectrum, SourceSpectrum
from synphot import units

from etc.utils import read_element, read_eso_spectra, percentage_difference


class TestReadElement:
    eso_unit = units.u.photon / units.u.s / units.u.m**2 / units.u.um

    def test_read_element(self):
        test_fp = files("tests.etc.data").joinpath("test_mirror.dat").as_posix()

        element = read_element(test_fp)

        assert_quantity_allclose(element.waveset[0], 2750 * u.AA)
        assert_quantity_allclose(element(element.waveset[0]), 0.64)

    def test_read_element_microns(self):
        test_fp = files("tests.etc.data").joinpath("test_mirror_microns.dat").as_posix()

        element = read_element(test_fp)

        assert_quantity_allclose(element.waveset[0], 3000 * u.AA)
        assert_quantity_allclose(element(element.waveset[0]), 0.751)

    def test_read_element_microns_passed_waveunit(self):
        test_fp = files("tests.etc.data").joinpath("test_mirror_microns.dat").as_posix()

        element = read_element(test_fp, wave_units=u.micron)

        assert_quantity_allclose(element.waveset[0], 3000 * u.AA)
        assert_quantity_allclose(element(element.waveset[0]), 0.751)

    def test_read_element_microns_passed_fluxunit(self):
        test_fp = files("tests.etc.data").joinpath("test_mirror_microns.dat").as_posix()

        element = read_element(test_fp, flux_units=u.percent)

        assert_quantity_allclose(element.waveset[0], 3000 * u.AA)
        assert_quantity_allclose(element(element.waveset[0]), 0.00751)  # File is not in percent

    def test_read_element_angstroms(self):
        test_fp = files("tests.etc.data").joinpath("test_atmos_angstroms.dat").as_posix()

        element = read_element(test_fp)

        assert_quantity_allclose(element.waveset[0], 3200 * u.AA)
        assert_quantity_allclose(element(element.waveset[0]), 0.306, rtol=1e-3)

    def test_read_optic(self):
        test_fp = files("tests.etc.data").joinpath("test_optic_nm.csv").as_posix()

        element = read_element(test_fp)

        assert_quantity_allclose(element.waveset[0], 200 * u.nm)
        assert_quantity_allclose(element(element.waveset[0]), 0.000301121 * u.percent)
        assert_quantity_allclose(element.waveset[-1], 2200 * u.nm)
        assert_quantity_allclose(element(element.waveset[-1]), 85.98 * u.percent)

    def test_read_optic_passed_waveunit(self):
        test_fp = files("tests.etc.data").joinpath("test_optic_nm.csv").as_posix()

        element = read_element(test_fp, wave_units=u.nm)

        assert_quantity_allclose(element.waveset[0], 200 * u.nm)
        assert_quantity_allclose(element(element.waveset[0]), 0.000301121 * u.percent)
        assert_quantity_allclose(element.waveset[-1], 2200 * u.nm)
        assert_quantity_allclose(element(element.waveset[-1]), 85.98 * u.percent)

    def test_read_optic_microns(self):
        test_fp = files("tests.etc.data").joinpath("test_optic_microns.csv").as_posix()

        element = read_element(test_fp)

        assert_quantity_allclose(element.waveset[0], 6 * u.AA)
        assert_quantity_allclose(element(element.waveset[0]), 0.0)
        assert_quantity_allclose(element.waveset[1], 6.08 * u.AA)
        assert_quantity_allclose(element(element.waveset[1]), 20 * u.percent)
        assert_quantity_allclose(element.waveset[-1], 130 * u.AA)
        assert_quantity_allclose(element(element.waveset[0]), 0.0)

    def test_read_optic_microns_passed_waveunit(self):
        test_fp = files("tests.etc.data").joinpath("test_optic_microns.csv").as_posix()

        element = read_element(test_fp, wave_units=u.micron)

        assert_quantity_allclose(element.waveset[0], 0.6 * u.micron)
        assert_quantity_allclose(element(element.waveset[0]), 0.0)
        assert_quantity_allclose(element.waveset[1], 0.608 * u.micron)
        assert_quantity_allclose(element(element.waveset[1]), 20 * u.percent)
        assert_quantity_allclose(element.waveset[-1], 13 * u.micron)
        assert_quantity_allclose(element(element.waveset[0]), 0.0)

    def test_read_element_as_sourcespectrum(self):
        test_fp = files("tests.etc.data").joinpath("test_radiance.dat").as_posix()

        element = read_element(test_fp, element_type="spectrum")

        assert isinstance(element, SourceSpectrum)
        assert element.waveset[0].unit == u.AA
        assert_quantity_allclose(element.waveset[0], 300 * u.nm)
        assert_quantity_allclose(element(element.waveset[0]), 16.4805 * units.PHOTLAM, rtol=1e-3)

    def test_read_element_as_sourcespectrum_ownfluxunits(self):
        test_fp = files("tests.etc.data").joinpath("test_radiance.dat").as_posix()

        element = read_element(test_fp, element_type="spectrum", flux_units=self.eso_unit)

        assert isinstance(element, SourceSpectrum)
        assert element.waveset[0].unit == u.AA
        assert element(element.waveset[0]).unit == units.PHOTLAM
        assert_quantity_allclose(element.waveset[0], 300 * u.nm)
        assert_quantity_allclose(element(element.waveset[0]), 16.4805 * self.eso_unit, rtol=1e-3)

    def test_read_fits_as_radiance(self):
        test_fp = files("tests.etc.data").joinpath("test_skytable.fits").as_posix()

        element = read_element(test_fp, element_type="radiance")

        assert isinstance(element, SourceSpectrum)
        assert element.waveset[0].unit == u.AA
        assert element(element.waveset[0]).unit == units.PHOTLAM
        assert_quantity_allclose(element.waveset[0], 300 * u.nm)
        assert_quantity_allclose(element.waveset[-1], 1200 * u.nm)
        assert_quantity_allclose(element(element.waveset[0]), 14.20246 * self.eso_unit, rtol=1e-3)
        assert_quantity_allclose(element(element.waveset[-1]), 3499.22056 * self.eso_unit, rtol=1e-3)

    def test_read_fits_as_radiance_ownfluxunits(self):
        test_fp = files("tests.etc.data").joinpath("test_skytable.fits").as_posix()

        element = read_element(test_fp, element_type="RADiANCE", flux_units=self.eso_unit)

        assert isinstance(element, SourceSpectrum)
        assert element.waveset[0].unit == u.AA
        assert element(element.waveset[0]).unit == units.PHOTLAM
        assert_quantity_allclose(element.waveset[0], 300 * u.nm)
        assert_quantity_allclose(element.waveset[-1], 1200 * u.nm)
        assert_quantity_allclose(element(element.waveset[0]), 14.20246 * self.eso_unit, rtol=1e-3)
        assert_quantity_allclose(element(element.waveset[-1]), 3499.22056 * self.eso_unit, rtol=1e-3)

    def test_read_fits_as_element(self):
        test_fp = files("tests.etc.data").joinpath("test_skytable.fits").as_posix()

        element = read_element(test_fp, element_type="spectral_element")

        assert isinstance(element, SpectralElement)
        assert element.waveset[0].unit == u.AA
        assert element(element.waveset[0]).unit == units.THROUGHPUT
        assert_quantity_allclose(element.waveset[0], 300 * u.nm)
        assert_quantity_allclose(element.waveset[-1], 1200 * u.nm)
        assert_quantity_allclose(element(element.waveset[0]), 0.03626, rtol=1e-3)
        assert_quantity_allclose(element(element.waveset[-1]), 0.971, rtol=1e-3)

    def test_read_hst_fits_as_element(self):
        test_fp = files("tests.etc.data").joinpath("test_ota.fits").as_posix()

        element = read_element(test_fp, element_type="spectral_element")

        assert isinstance(element, SpectralElement)
        assert element.waveset[0].unit == u.AA
        assert element(element.waveset[0]).unit == units.THROUGHPUT
        assert_quantity_allclose(element.waveset[0], 500 * u.AA)
        assert_quantity_allclose(element.waveset[1], 1000 * u.AA)
        assert_quantity_allclose(element.waveset[-2], 30000 * u.AA)
        assert_quantity_allclose(element(element.waveset[0]), 0.0, rtol=1e-3)
        assert_quantity_allclose(element(element.waveset[1]), 0.023, rtol=1e-3)
        assert_quantity_allclose(element(element.waveset[-2]), 0.731268, rtol=1e-3)

    def test_read_fudge_greater_than_one(self):
        """ESO Omegacam has an optics throughput fudge which goes fro, 0.93 to 3.2
        with a mean of 1.4024. This incorrectly triggers the divide by 100 logic
        to convert it from a percentage to a 0..1 range. This test checks this
        is correctly handled"""
        test_fp = files("tests.etc.data").joinpath("test_fudge.dat").as_posix()

        element = read_element(test_fp)

        assert_quantity_allclose(element.waveset[0], 320 * u.nm)
        assert_quantity_allclose(element(element.waveset[0]), 0.93)
        assert_quantity_allclose(element.waveset[-1], 1000 * u.nm)
        assert_quantity_allclose(element(element.waveset[-1]), 3.2)

    def test_read_fudge_greater_than_one_as_SE(self):
        test_fp = files("tests.etc.data").joinpath("test_fudge.dat").as_posix()

        element = read_element(test_fp, element_type="spectral_element")

        assert type(element) == SpectralElement
        assert_quantity_allclose(element.waveset[0], 320 * u.nm)
        assert_quantity_allclose(element(element.waveset[0]), 0.93)
        assert_quantity_allclose(element.waveset[-1], 1000 * u.nm)
        assert_quantity_allclose(element(element.waveset[-1]), 3.2)

    def test_config_item_Vfilter(self):
        element = read_element("V")

        assert type(element) == SpectralElement
        assert_quantity_allclose(element.waveset[0], 4700 * u.AA)
        assert_quantity_allclose(element(element.waveset[0]), 0.0)
        assert_quantity_allclose(element.waveset[-1], 7000 * u.AA)
        assert_quantity_allclose(element(element.waveset[-1]), 0.0)
        assert_quantity_allclose(element(530 * u.nm), 1.0)


class TestReadESOSpectra:
    @classmethod
    def setup_class(cls):
        cls.test_fits_fp = files("tests.etc.data").joinpath("test_spectrum.fits").as_posix()
        cls.funit = u.erg / (u.cm**2 * u.s * u.AA)

    def test_spectra(self):
        source_spec = read_eso_spectra(self.test_fits_fp)

        assert type(source_spec) == SourceSpectrum
        assert_quantity_allclose(source_spec.waveset[0], 3451.65 * u.AA)
        assert_quantity_allclose(
            source_spec(source_spec.waveset[0], flux_unit="FLAM"), 4.341547e-16 * self.funit
        )


class TestPercentageDifference:
    def test1(self):
        expected_value = 0.817881 * u.percent

        result = percentage_difference(149.03496401806, 147.821)

        assert_quantity_allclose(expected_value, result, rtol=1e-6)

    def test2(self):
        expected_value = 0.817881 * u.percent

        result = percentage_difference(147.821, 149.03496401806)

        assert_quantity_allclose(expected_value, result, rtol=1e-6)
