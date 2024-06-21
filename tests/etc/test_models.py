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
from synphot import units
from synphot.spectrum import SpectralElement, BaseUnitlessSpectrum, SourceSpectrum

from etc.models import *
from etc.utils import ETCError


class TestSite:
    @pytest.fixture(autouse=True)
    def setUp(self):
        self.eso_rad_unit = u.photon / u.s / u.m**2 / u.um

        self.zp_B = (4063 * u.Jy).to(units.PHOTLAM, equivalencies=u.spectral_density(4300 * u.AA))
        self.zp_V = (3636 * u.Jy).to(units.PHOTLAM, equivalencies=u.spectral_density(5500 * u.AA))

    def test_initialize_defaults(self):
        site = Site()

        assert site.name == "Undefined"
        assert site.altitude == None
        assert site.latitude == None
        assert site.longitude == None

    def test_initialize1(self):
        site = Site(name="Maui", altitude=3055.0, latitude=20.7, longitude=-156.2)

        assert site.name == "Maui"
        assert site.altitude == 3055 * u.m
        assert site.latitude == 20.7 * u.deg
        assert site.longitude == -156.2 * u.deg

    def test_initialize2(self):
        test_config_file = toml.load(files("tests.etc.data").joinpath("test1.toml"))
        site = Site(**test_config_file["site"])

        assert site.name == "BPL"
        assert site.altitude == 17 * u.m
        assert site.latitude == 34.5 * u.deg
        assert site.longitude == -119.86 * u.deg
        assert site.tpeak() == 0.9

    def test_transmission_file(self):
        test_config = {
            "name": "BPL",
            "altitude": 7,
            "latitude": 34.5,
            "longitude": -119.86,
            "transmission": files("tests.etc.data").joinpath("test_atmos.dat").as_posix(),
        }
        site = Site(**test_config)

        assert site.name == "BPL"
        assert site.altitude == 7 * u.m
        assert site.latitude == 34.5 * u.deg
        assert site.longitude == -119.86 * u.deg
        assert site.tpeak() == 0.975

    def test_radiance_file(self):
        test_config = {
            "name": "Cerro Armazones",
            "altitude": 3060,
            "latitude": -24.58916667,
            "longitude": -70.19166667,
            "radiance": files("tests.etc.data").joinpath("test_radiance.dat").as_posix(),
        }
        site = Site(**test_config)

        assert site.name == "Cerro Armazones"
        assert site.altitude == 3060 * u.m
        assert site.latitude == -24.58916667 * u.deg
        assert site.longitude == -70.19166667 * u.deg
        radiance = site.radiance
        assert isinstance(radiance, SourceSpectrum)
        assert_quantity_allclose(radiance.waveset[1], 350 * u.nm)
        assert_quantity_allclose(radiance(radiance.waveset[1]), 371.148 * self.eso_rad_unit)

    def test_radiance_file_bad_units(self):
        test_config = {
            "name": "Cerro Armazones",
            "altitude": 3060,
            "latitude": -24.58916667,
            "longitude": -70.19166667,
            "radiance": files("tests.etc.data").joinpath("test_radiance.dat").as_posix(),
            "radiance_units": "photon/elephant**2",
        }
        site = Site(**test_config)

        assert site.name == "Cerro Armazones"
        assert site.altitude == 3060 * u.m
        assert site.latitude == -24.58916667 * u.deg
        assert site.longitude == -70.19166667 * u.deg
        radiance = site.radiance
        assert isinstance(radiance, SourceSpectrum)
        assert_quantity_allclose(radiance.waveset[1], 350 * u.nm)
        assert_quantity_allclose(radiance(radiance.waveset[1]), 371.148 * self.eso_rad_unit)

    def test_radiance_file_goodunits(self):
        test_config = {
            "name": "Cerro Armazones",
            "altitude": 3060,
            "latitude": -24.58916667,
            "longitude": -70.19166667,
            "radiance": files("tests.etc.data").joinpath("test_radiance.dat").as_posix(),
            "radiance_units": "photon/s/m**2/nm",
        }
        site = Site(**test_config)

        assert site.name == "Cerro Armazones"
        assert site.altitude == 3060 * u.m
        assert site.latitude == -24.58916667 * u.deg
        assert site.longitude == -70.19166667 * u.deg
        radiance = site.radiance
        assert isinstance(radiance, SourceSpectrum)
        assert_quantity_allclose(radiance.waveset[1], 350 * u.nm)
        # 1000x bigger due to micron->nm
        assert_quantity_allclose(radiance(radiance.waveset[1]), 371148 * self.eso_rad_unit)

    def test_sky_spectrum_default(self):
        test_config = {
            "name": "Cerro Armazones",
            "altitude": 3060,
            "latitude": -24.58916667,
            "longitude": -70.19166667,
        }
        site = Site(**test_config)

        sky_spectrum = site.sky_spectrum()
        assert isinstance(sky_spectrum, SourceSpectrum)
        assert_quantity_allclose(sky_spectrum.waveset[0], 300 * u.nm)
        assert_quantity_allclose(sky_spectrum.waveset[-1], 1200 * u.nm)
        assert_quantity_allclose(sky_spectrum(sky_spectrum.waveset[0]), self.zp_V)

    def test_sky_spectrum_V(self):
        test_config = {
            "name": "Cerro Armazones",
            "altitude": 3060,
            "latitude": -24.58916667,
            "longitude": -70.19166667,
        }
        site = Site(**test_config)

        sky_spectrum = site.sky_spectrum("V")
        assert isinstance(sky_spectrum, SourceSpectrum)
        assert_quantity_allclose(sky_spectrum.waveset[0], 300 * u.nm)
        assert_quantity_allclose(sky_spectrum.waveset[-1], 1200 * u.nm)
        assert_quantity_allclose(sky_spectrum(sky_spectrum.waveset[0]), self.zp_V)

    def test_sky_spectrum_B(self):
        test_config = {
            "name": "Cerro Armazones",
            "altitude": 3060,
            "latitude": -24.58916667,
            "longitude": -70.19166667,
        }
        site = Site(**test_config)

        sky_spectrum = site.sky_spectrum("B")
        assert isinstance(sky_spectrum, SourceSpectrum)
        assert_quantity_allclose(sky_spectrum.waveset[0], 300 * u.nm)
        assert_quantity_allclose(sky_spectrum.waveset[-1], 1200 * u.nm)
        assert_quantity_allclose(sky_spectrum(sky_spectrum.waveset[0]), self.zp_B)

    def test_sky_spectrum_radiance_file(self):
        test_config = {
            "name": "Cerro Armazones",
            "altitude": 3060,
            "latitude": -24.58916667,
            "longitude": -70.19166667,
            "radiance": files("tests.etc.data").joinpath("test_radiance.dat").as_posix(),
        }
        site = Site(**test_config)

        sky_spectrum = site.sky_spectrum()
        assert isinstance(sky_spectrum, SourceSpectrum)
        assert_quantity_allclose(sky_spectrum.waveset[1], 350 * u.nm)
        assert_quantity_allclose(sky_spectrum(sky_spectrum.waveset[1]), 371.148 * self.eso_rad_unit)

    def test_rebin_transmission(self):
        test_config = {
            "name": "BPL",
            "altitude": 7,
            "latitude": 34.5,
            "longitude": -119.86,
            "transmission": files("tests.etc.data").joinpath("test_atmos.dat").as_posix(),
        }
        site = Site(**test_config)
        new_waves = u.Quantity([300, 600, 900, 1200], u.nm)

        new_trans = site.rebin_transmission(new_waves)
        assert isinstance(new_trans, SpectralElement)

    def test_sky_brightness(self):
        expected_sky_mags = {"Days_from_New_Moon": 0, "U": 22.0, "B": 22.7, "V": 21.8, "R": 20.9, "I": 19.9}

        site = Site()

        assert expected_sky_mags == site.sky_mags

    def test_extinction_to_transmission(self):
        expected_transmission = 0.9

        site = Site()
        extinction = 0.114393726402

        transmission = site._extinction_to_transmission(extinction)

        assert_quantity_allclose(expected_transmission, transmission)

    def test_extinction_to_transmission_B_a1p5(self):
        expected_transmission = 0.707945764

        site = Site()
        extinction = 0.25
        airmass = 1.5

        transmission = site._extinction_to_transmission(extinction, airmass)

        assert_quantity_allclose(expected_transmission, transmission)

    def test_extinction_to_transmission_Z_a2p0(self):
        expected_transmission = 0.9120108393559098

        site = Site()
        extinction = 0.05
        airmass = 2.0

        transmission = site._extinction_to_transmission(extinction, airmass)

        assert_quantity_allclose(expected_transmission, transmission)

    def test_transmission_to_extinction(self):
        expected_extinct = 0.114393726402

        site = Site()
        transmission = 0.9

        extinct = site._transmission_to_extinction(transmission)

        assert_quantity_allclose(expected_extinct, extinct)

    def test_transmission_to_extinction_B_a1p5(self):
        expected_extinct = 0.25

        site = Site()
        transmission = 0.707945764
        airmass = 1.5

        extinct = site._transmission_to_extinction(transmission, airmass)

        assert_quantity_allclose(expected_extinct, extinct)

    def test_transmission_to_extinction_Z_a2p0(self):
        expected_extinct = 0.05

        site = Site()
        transmission = 0.9120108393559098
        airmass = 2.0

        extinct = site._transmission_to_extinction(transmission, airmass)

        assert_quantity_allclose(expected_extinct, extinct)


class TestTelescope:
    def test_initialize_defaults(self):
        tel = Telescope()

        assert tel.name == "Undefined"
        assert tel.size == 0 * u.m
        assert tel.area == 0 * u.m * u.m
        assert tel.num_mirrors == 2

    def test_initialize1(self):
        tel = Telescope(name="FTN", size=2.0, area=2.574, num_mirrors=3)

        assert tel.name == "FTN"
        assert tel.size == 2 * u.m
        assert tel.area == 2.574 * u.m * u.m
        assert tel.num_mirrors == 3
        assert tel.tpeak() == 0.85**3

    def test_initialize2(self):
        test_config_file = toml.load(files("tests.etc.data").joinpath("test1.toml"))
        tel = Telescope(**test_config_file["telescope"])

        assert tel.name == "BPL 1-m"
        assert tel.size == 1 * u.m
        assert tel.area == 0.625 * u.m * u.m
        assert tel.num_mirrors == 2
        assert tel.tpeak() == 0.8 * 0.8

    def test_reflectivity_file(self):
        test_config = {
            "name": "BPL 1-m",
            "size": 1,
            "area": 0.625,
            "reflectivity": files("tests.etc.data").joinpath("test_mirror.dat").as_posix(),
        }
        tel = Telescope(**test_config)

        assert tel.name == "BPL 1-m"
        assert tel.size == 1 * u.m
        assert tel.area == 0.625 * u.m * u.m
        assert tel.num_mirrors == 2
        assert tel.tpeak() == 0.92**2

    def test_reflectivity_list(self):
        test_fp = files("tests.etc.data").joinpath("test_mirror.dat").as_posix()

        test_config = {"name": "BPL 1-m", "size": 1, "area": 0.625, "reflectivity": [test_fp, test_fp]}
        tel = Telescope(**test_config)

        assert tel.name == "BPL 1-m"
        assert tel.size == 1 * u.m
        assert tel.area == 0.625 * u.m * u.m
        assert tel.num_mirrors == 2
        assert tel.tpeak() == 0.92**2


class TestInstrument:
    def test_initialize_defaults(self):
        inst = Instrument()

        assert inst.name == "Undefined"
        assert inst.inst_type == "IMAGER"

    def test_initialize_bad_insttype(self):
        inst = Instrument(inst_type="Wibble")

        assert inst.name == "Undefined"
        assert inst.inst_type == "IMAGER"

    def test_initialize_spectrograph(self):
        inst = Instrument(name="FLOYDS", inst_type="Spectrograph")

        assert inst.name == "FLOYDS"
        assert inst.inst_type == "SPECTROGRAPH"
        assert inst.is_fiberfed == False

    def test_initialize_fiber_spectrograph(self):
        inst = Instrument(name="FEROS", inst_type="Spectrograph", fiber_diameter=2 * u.arcsec)

        assert inst.name == "FEROS"
        assert inst.inst_type == "SPECTROGRAPH"
        assert inst.is_fiberfed == True

    def test_trans_default(self):
        inst = Instrument()

        assert inst.num_mirrors == 0
        assert inst.num_lenses == 1
        assert inst.num_ar_coatings == 2
        assert inst.filterlist == []

        assert inst.transmission.tpeak() == 0.911493  # 0.93 (lens) * 0.99^2 (AR)

    def test_trans_modify_optics(self):
        optics_options = {"inst_lens_trans": 0.85, "inst_ar_coating_refl": 0.95}
        inst = Instrument(**optics_options)

        assert inst.num_mirrors == 0
        assert inst.num_lenses == 1
        assert inst.num_ar_coatings == 2

        assert_quantity_allclose(inst.transmission.tpeak(), 0.767125, 1e-5)

    def test_trans_spectrograph(self):
        optics_options = {
            "num_ar_coatings": 6,  # Air-glass interfaces: prism (2 sides), field flattener (4 sides)
            "num_inst_lenses": 3,  # Fused silica prism (two passes) and CCD window
            "num_inst_mirrors": 4,  # Mirrors:  Collimator, Fold, Camera, Grating
        }
        inst = Instrument(name="FLOYDS", inst_type="SPECTROGRAPH", **optics_options)

        assert inst.name == "FLOYDS"
        assert inst.inst_type == "SPECTROGRAPH"
        assert inst.num_mirrors == 4
        assert inst.num_lenses == 3
        assert inst.num_ar_coatings == 6

        assert_quantity_allclose(inst.transmission.tpeak(), 0.73482187, 1e-5)

    def test_filterset(self):
        optics_options = {
            "filterlist": ["u", "g", "r", "i", "z"],
        }
        inst = Instrument(**optics_options)

        assert inst.filterlist == ["u", "g", "r", "i", "z"]
        assert len(inst.filterset) == 5
        assert isinstance(inst.filterset["g"], SpectralElement)

    def test_filterset_primenotation(self):
        optics_options = {
            "filterlist": ["up", "gp", "rp", "ip", "zs"],
        }
        inst = Instrument(**optics_options)

        assert inst.filterlist == ["up", "gp", "rp", "ip", "zs"]
        assert len(inst.filterset) == 5
        assert isinstance(inst.filterset["gp"], SpectralElement)
        assert isinstance(inst.filterset["zs"], SpectralElement)

    def test_filterset_mixedcase(self):
        optics_options = {
            "filterlist": ["r", "R"],
        }

        inst = Instrument(**optics_options)

        assert inst.filterlist == ["r", "R"]
        assert len(inst.filterset) == 2
        assert isinstance(inst.filterset["r"], SpectralElement)
        assert isinstance(inst.filterset["R"], SpectralElement)
        assert inst.filterset["r"].meta["header"]["descrip"] != inst.filterset["R"].meta["header"]["descrip"]
        assert inst.filterset["r"].meta["expr"] != inst.filterset["R"].meta["expr"]

    def test_throughput(self):
        optics_options = {
            "filterlist": [
                "r",
            ]
        }

        inst = Instrument(**optics_options)

        assert inst.filterlist == [
            "r",
        ]
        assert len(inst.filterset) == 1
        assert isinstance(inst.filterset["r"], SpectralElement)
        assert isinstance(inst.throughput("r"), SpectralElement)
        # 0.93 (lens) * 0.99^2 (AR) * 0.99165802 (filter) * 0.9 (CCD)
        assert_quantity_allclose(inst.throughput("r").tpeak(), 0.81350041, 1e-5)

    def test_throughput_new_ccd_qe(self):
        optics_options = {
            "filterlist": [
                "r",
            ],
            "ccd_qe": 0.45,
        }

        inst = Instrument(**optics_options)

        assert inst.filterlist == [
            "r",
        ]
        assert len(inst.filterset) == 1
        assert isinstance(inst.filterset["r"], SpectralElement)
        assert isinstance(inst.throughput("r"), SpectralElement)
        assert_quantity_allclose(inst.throughput("r").tpeak(), 0.81350041 / 2.0, 1e-5)

    def test_throughput_ccd_qe_quantity(self):
        optics_options = {
            "filterlist": [
                "r",
            ],
            "ccd_qe": u.Quantity(0.45, u.dimensionless_unscaled),
        }

        inst = Instrument(**optics_options)

        assert inst.filterlist == [
            "r",
        ]
        assert len(inst.filterset) == 1
        assert isinstance(inst.filterset["r"], SpectralElement)
        assert isinstance(inst.throughput("r"), SpectralElement)
        assert_quantity_allclose(inst.throughput("r").tpeak(), 0.81350041 / 2.0, 1e-5)

    def test_throughput_ccd_qe_file(self):
        optics_options = {
            "filterlist": [
                "r",
            ],
            "ccd_qe": files("tests.etc.data").joinpath("test_ccd_qe.dat"),
        }

        inst = Instrument(**optics_options)

        assert inst.filterlist == [
            "r",
        ]
        assert len(inst.filterset) == 1
        assert isinstance(inst.filterset["r"], SpectralElement)
        assert isinstance(inst.throughput("r"), SpectralElement)
        assert isinstance(inst.ccd_qe, BaseUnitlessSpectrum)
        assert_quantity_allclose(inst.throughput("r").tpeak(), 0.8856121333333334, 1e-5)

    def test_throughput_ccd_qe_dict(self):
        optics_options = {
            "filterlist": ["U", "B", "V", "R", "I"],
            "ccd_qe": {"U": 0.24, "B": 0.44, "V": 0.7, "R": 0.79, "I": 0.59},
        }

        inst = Instrument(**optics_options)

        assert len(inst.filterset) == 5
        assert inst.filterlist == ["U", "B", "V", "R", "I"]
        assert isinstance(inst.filterset["U"], SpectralElement)
        assert isinstance(inst.throughput("U"), SpectralElement)
        assert_quantity_allclose(
            inst.throughput("U").tpeak(), 1 * 0.93 * 0.99**2 * optics_options["ccd_qe"]["U"], 1e-5
        )

    def test_ccd_parameters(self):
        expected_gain = 0.9 * (u.electron / u.adu)
        expected_noise = 4.0 * (u.electron / u.pix)

        ccd_options = {
            "ccd_readnoise": 4.0,
            "ccd_gain": 0.9,
            "ccd_qe": files("tests.etc.data").joinpath("test_ccd_qe.dat"),
        }

        inst = Instrument(**ccd_options)

        assert isinstance(inst.ccd_qe, BaseUnitlessSpectrum)
        assert_quantity_allclose(inst.ccd_gain, expected_gain)
        assert_quantity_allclose(inst.ccd_readnoise, expected_noise)

    def test_ccdpixsize_default(self):
        expected_value = 17.5 * u.micron

        inst_options = {"ccd_pixsize": 17.5}

        inst = Instrument(**inst_options)

        assert isinstance(inst.ccd_pixsize, u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_pixsize)

    def test_ccdpixsize_goodunits(self):
        expected_value = 13.5 * u.micron

        inst_options = {"ccd_pixsize": 0.0135, "ccd_pixsize_units": "mm"}

        inst = Instrument(**inst_options)

        assert isinstance(inst.ccd_pixsize, u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_pixsize)

    def test_ccdpixsize_goodunits2(self):
        expected_value = 13.5 * u.micron

        inst_options = {"ccd_pixsize": 13.5e-6, "ccd_pixsize_units": "m"}

        inst = Instrument(**inst_options)

        assert isinstance(inst.ccd_pixsize, u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_pixsize)

    def test_ccdpixsize_badunits(self):
        expected_value = 13.5 * u.micron

        inst_options = {"ccd_pixsize": 13.5, "ccd_pixsize_units": "wibble"}

        inst = Instrument(**inst_options)

        assert isinstance(inst.ccd_pixsize, u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_pixsize)

    def test_focalscale_default(self):
        expected_value = 17.5 * (u.arcsec / u.mm)

        inst_options = {"focal_scale": 17.5}

        inst = Instrument(**inst_options)

        assert isinstance(inst.focal_scale, u.Quantity)
        assert_quantity_allclose(expected_value, inst.focal_scale)

    def test_focalscale_goodunits(self):
        expected_value = 1.75 * (u.arcsec / u.cm)

        inst_options = {"focal_scale": 1.75, "focal_scale_units": "arcsec/cm"}

        inst = Instrument(**inst_options)

        assert isinstance(inst.focal_scale, u.Quantity)
        assert_quantity_allclose(expected_value, inst.focal_scale)

    def test_focalscale_badunits(self):
        expected_value = 17.5 * (u.arcsec / u.mm)

        inst_options = {"focal_scale": 17.5, "focal_scale_units": "wibble/cm"}

        inst = Instrument(**inst_options)

        assert isinstance(inst.focal_scale, u.Quantity)
        assert_quantity_allclose(expected_value, inst.focal_scale)

    def test_ccdpixscale(self):
        expected_value = 0.23625 * u.arcsec

        inst_options = {"focal_scale": 17.5, "ccd_pixsize": 13.5}

        inst = Instrument(**inst_options)

        assert isinstance(inst.focal_scale, u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_pixscale)

    def test_ccdpixscale_goodunits(self):
        expected_value = 0.23625 * u.arcsec

        inst_options = {
            "focal_scale": 175,
            "focal_scale_units": "arcsec/cm",
            "ccd_pixsize": 0.0135,
            "ccd_pixsize_units": "mm",
        }

        inst = Instrument(**inst_options)

        assert isinstance(inst.focal_scale, u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_pixscale)

    def test_ccdpixscale_badunits(self):
        expected_value = 0.23625 * u.arcsec

        inst_options = {
            "ccd_pixsize": 13.5,
            "ccd_pixsize_units": "wibble",
            "focal_scale": 17.5,
            "focal_scale_units": "more/wibble",
        }

        inst = Instrument(**inst_options)

        assert isinstance(inst.focal_scale, u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_pixscale)

    def test_ccdpixscale_binning(self):
        expected_value = 0.23625 * 2 * u.arcsec

        inst_options = {"focal_scale": 17.5, "ccd_pixsize": 13.5, "ccd_xbinning": 2, "ccd_ybinning": 2}

        inst = Instrument(**inst_options)

        assert isinstance(inst.focal_scale, u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_pixscale)

    def test_ccdpixscale_badbinning(self):
        inst_options = {"focal_scale": 17.5, "ccd_pixsize": 13.5, "ccd_xbinning": 2, "ccd_ybinning": 1}

        with pytest.raises(Exception) as execinfo:
            inst = Instrument(**inst_options)

    def test_ccdfov_square_defaults(self):
        expected_value = (483.84 * u.arcsec, 483.84 * u.arcsec)

        inst_options = {
            "ccd_pixsize": 13.5,
            "focal_scale": 17.5,
            "ccd_xpixels": 2048,
            "ccd_ypixels": 2048,
        }

        inst = Instrument(**inst_options)

        assert isinstance(inst.ccd_fov()[0], u.Quantity)
        assert isinstance(inst.ccd_fov()[1], u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_fov())

    def test_ccdfov_square_arcmin(self):
        expected_value = (483.84 / 60.0 * u.arcmin, 483.84 / 60.0 * u.arcmin)

        inst_options = {
            "ccd_pixsize": 13.5,
            "focal_scale": 17.5,
            "ccd_xpixels": 2048,
            "ccd_ypixels": 2048,
        }

        inst = Instrument(**inst_options)

        fov = inst.ccd_fov(u.arcmin)
        assert isinstance(fov[0], u.Quantity)
        assert isinstance(fov[1], u.Quantity)
        assert_quantity_allclose(expected_value, fov)

    def test_ccdfov_square_badunits_output(self):
        inst_options = {
            "ccd_pixsize": 13.5,
            "focal_scale": 17.5,
            "ccd_xpixels": 2048,
            "ccd_ypixels": 2048,
        }

        inst = Instrument(**inst_options)

        with pytest.raises(Exception) as execinfo:
            fov = inst.ccd_fov(u.m)

    def test_ccdfov_square_mm(self):
        expected_value = (483.84 * u.arcsec, 483.84 * u.arcsec)

        inst_options = {
            "ccd_pixsize": 0.0135,
            "ccd_pixsize_units": "mm",
            "focal_scale": 17.5,
            "ccd_xpixels": 2048,
            "ccd_ypixels": 2048,
        }

        inst = Instrument(**inst_options)

        assert isinstance(inst.ccd_fov()[0], u.Quantity)
        assert isinstance(inst.ccd_fov()[1], u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_fov())

    def test_ccdfov_square_mm_cm(self):
        expected_value = (483.84 * u.arcsec, 483.84 * u.arcsec)

        inst_options = {
            "ccd_pixsize": 0.0135,
            "ccd_pixsize_units": "mm",
            "focal_scale": 175,
            "focal_scale_units": "arcsec/cm",
            "ccd_xpixels": 2048,
            "ccd_ypixels": 2048,
        }

        inst = Instrument(**inst_options)

        assert isinstance(inst.ccd_fov()[0], u.Quantity)
        assert isinstance(inst.ccd_fov()[1], u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_fov())

    def test_ccdfov_rect_defaults(self):
        expected_value = (1754 * u.arcsec, 1169.3 * u.arcsec)

        inst_options = {
            "ccd_pixsize": 9,
            "ccd_pixsize_units": "um",
            "focal_scale": 63.44,
            "focal_scale_units": "arcsec/mm",
            "ccd_xpixels": 3072,
            "ccd_ypixels": 2048,
        }

        inst = Instrument(**inst_options)

        assert isinstance(inst.ccd_fov()[0], u.Quantity)
        assert isinstance(inst.ccd_fov()[1], u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_fov(), rtol=1e-1)

    def test_ccdfov_rect_arcmin(self):
        expected_value = (29.2 * u.arcmin, 19.5 * u.arcmin)

        inst_options = {
            "ccd_pixsize": 9,
            "ccd_pixsize_units": "um",
            "focal_scale": 63.44,
            "focal_scale_units": "arcsec/mm",
            "ccd_xpixels": 3072,
            "ccd_ypixels": 2048,
        }

        inst = Instrument(**inst_options)

        assert isinstance(inst.ccd_fov()[0], u.Quantity)
        assert isinstance(inst.ccd_fov()[1], u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_fov(u.arcmin), rtol=1e-1)

    def test_ccdfov_max_fov(self):
        size = (7 * u.arcmin).to(u.arcsec)
        expected_value = (size, size)

        inst_options = {
            "ccd_pixsize": 10,
            "ccd_pixsize_units": "um",
            "focal_scale": 10.12,
            "focal_scale_units": "arcsec/mm",
            "ccd_xpixels": 8120,
            "ccd_ypixels": 8120,
            "fov_xsize": 7,
            "fov_ysize": 7,
        }

        inst = Instrument(**inst_options)

        assert isinstance(inst.ccd_fov()[0], u.Quantity)
        assert isinstance(inst.ccd_fov()[1], u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_fov(), rtol=1e-1)

    def test_ccdfov_max_fov_arcmin(self):
        size = 7 * u.arcmin
        expected_value = (size, size)

        inst_options = {
            "ccd_pixsize": 10,
            "ccd_pixsize_units": "um",
            "focal_scale": 10.12,
            "focal_scale_units": "arcsec/mm",
            "ccd_xpixels": 8120,
            "ccd_ypixels": 8120,
            "fov_xsize": 7,
            "fov_ysize": 7,
        }

        inst = Instrument(**inst_options)

        assert isinstance(inst.ccd_fov()[0], u.Quantity)
        assert isinstance(inst.ccd_fov()[1], u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_fov(u.arcmin), rtol=1e-1)

    def test_ccdfov_max_fov_units_arcmin(self):
        size = 7 * u.arcmin
        expected_value = (size, size)

        inst_options = {
            "ccd_pixsize": 10,
            "ccd_pixsize_units": "um",
            "focal_scale": 10.12,
            "focal_scale_units": "arcsec/mm",
            "ccd_xpixels": 8120,
            "ccd_ypixels": 8120,
            "fov_xsize": 420,
            "fov_ysize": 420,
            "fov_units": "arcsec",
        }

        inst = Instrument(**inst_options)

        assert isinstance(inst.ccd_fov()[0], u.Quantity)
        assert isinstance(inst.ccd_fov()[1], u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_fov(u.arcmin), rtol=1e-1)

    def test_ccdfov_max_fov_units_bad(self):
        size = 7 * u.arcmin
        expected_value = (size, size)

        inst_options = {
            "ccd_pixsize": 10,
            "ccd_pixsize_units": "um",
            "focal_scale": 10.12,
            "focal_scale_units": "arcsec/mm",
            "ccd_xpixels": 8120,
            "ccd_ypixels": 8120,
            "fov_xsize": 7,
            "fov_ysize": 7,
            "fov_units": "furlongs",
        }

        inst = Instrument(**inst_options)

        assert isinstance(inst.ccd_fov()[0], u.Quantity)
        assert isinstance(inst.ccd_fov()[1], u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_fov(u.arcmin), rtol=1e-1)

    def test_ccdbinning_defaults(self):
        expected_value = (1, 1)

        inst_options = {
            "ccd_pixsize": 10,
            "ccd_pixsize_units": "um",
            "focal_scale": 10.12,
            "focal_scale_units": "arcsec/mm",
            "ccd_xpixels": 8120,
            "ccd_ypixels": 8120,
        }

        inst = Instrument(**inst_options)

        assert expected_value == inst.ccd_binning

    def test_ccdbinning_bin1(self):
        expected_value = (1, 1)

        inst_options = {
            "ccd_pixsize": 10,
            "ccd_pixsize_units": "um",
            "focal_scale": 10.12,
            "focal_scale_units": "arcsec/mm",
            "ccd_xpixels": 8120,
            "ccd_ypixels": 8120,
            "ccd_xbinning": 1,
            "ccd_ybinning": 1,
        }

        inst = Instrument(**inst_options)

        assert expected_value == inst.ccd_binning

    def test_ccdbinning_bin2(self):
        expected_value = (2, 2)

        inst_options = {
            "ccd_pixsize": 10,
            "ccd_pixsize_units": "um",
            "focal_scale": 10.12,
            "focal_scale_units": "arcsec/mm",
            "ccd_xpixels": 8120,
            "ccd_ypixels": 8120,
            "ccd_xbinning": 2,
            "ccd_ybinning": 2,
        }

        inst = Instrument(**inst_options)

        assert expected_value == inst.ccd_binning

    def test_ccdbinning_bin3(self):
        expected_value = (3, 3)

        inst_options = {
            "ccd_pixsize": 10,
            "ccd_pixsize_units": "um",
            "focal_scale": 10.12,
            "focal_scale_units": "arcsec/mm",
            "ccd_xpixels": 8120,
            "ccd_ypixels": 8120,
            "ccd_xbinning": 3,
            "ccd_ybinning": 3,
        }

        inst = Instrument(**inst_options)

        assert expected_value == inst.ccd_binning

    def test_ccdbinning_badvalue_neg(self):
        expected_value = (1, 1)

        inst_options = {
            "ccd_pixsize": 10,
            "ccd_pixsize_units": "um",
            "focal_scale": 10.12,
            "focal_scale_units": "arcsec/mm",
            "ccd_xpixels": 8120,
            "ccd_ypixels": 8120,
            "ccd_xbinning": -3,
            "ccd_ybinning": -3,
        }

        inst = Instrument(**inst_options)

        assert expected_value == inst.ccd_binning

    def test_ccdbinning_badvalue_nonint(self):
        inst_options = {
            "ccd_pixsize": 10,
            "ccd_pixsize_units": "um",
            "focal_scale": 10.12,
            "focal_scale_units": "arcsec/mm",
            "ccd_xpixels": 8120,
            "ccd_ypixels": 8120,
            "ccd_xbinning": "foo",
            "ccd_ybinning": "bar",
        }

        with pytest.raises(ETCError) as execinfo:
            inst = Instrument(**inst_options)

    def test_ccdbinning_badvalue_nonint_y(self):
        inst_options = {
            "ccd_pixsize": 10,
            "ccd_pixsize_units": "um",
            "focal_scale": 10.12,
            "focal_scale_units": "arcsec/mm",
            "ccd_xpixels": 8120,
            "ccd_ypixels": 8120,
            "ccd_xbinning": 1,
            "ccd_ybinning": "bar",
        }

        with pytest.raises(ETCError) as execinfo:
            inst = Instrument(**inst_options)

    def test_ccdnumpixels_defaults(self):
        expected_value = (8120, 8120)

        inst_options = {
            "ccd_pixsize": 10,
            "ccd_pixsize_units": "um",
            "focal_scale": 10.12,
            "focal_scale_units": "arcsec/mm",
            "ccd_xpixels": 8120,
            "ccd_ypixels": 8120,
        }

        inst = Instrument(**inst_options)

        assert expected_value == inst.ccd_numpixels

    def test_ccdnumpixels_bin1(self):
        expected_value = (8120, 8120)

        inst_options = {
            "ccd_pixsize": 10,
            "ccd_pixsize_units": "um",
            "focal_scale": 10.12,
            "focal_scale_units": "arcsec/mm",
            "ccd_xpixels": 8120,
            "ccd_ypixels": 8120,
            "ccd_xbinning": 1,
            "ccd_ybinning": 1,
        }

        inst = Instrument(**inst_options)

        assert expected_value == inst.ccd_numpixels

    def test_ccdnumpixels_bin2(self):
        expected_value = (4060, 4060)

        inst_options = {
            "ccd_pixsize": 10,
            "ccd_pixsize_units": "um",
            "focal_scale": 10.12,
            "focal_scale_units": "arcsec/mm",
            "ccd_xpixels": 8120,
            "ccd_ypixels": 8120,
            "ccd_xbinning": 2,
            "ccd_ybinning": 2,
        }

        inst = Instrument(**inst_options)

        assert expected_value == inst.ccd_numpixels

    def test_ccdnumpixels_bin3(self):
        expected_value = (2706, 2706)

        inst_options = {
            "ccd_pixsize": 10,
            "ccd_pixsize_units": "um",
            "focal_scale": 10.12,
            "focal_scale_units": "arcsec/mm",
            "ccd_xpixels": 8120,
            "ccd_ypixels": 8120,
            "ccd_xbinning": 3,
            "ccd_ybinning": 3,
        }

        inst = Instrument(**inst_options)

        assert expected_value == inst.ccd_numpixels

    def test_vignette_imager(self):
        expected_vign = 1.0

        inst_options = {
            "inst_type": "imager",
        }

        inst = Instrument(**inst_options)

        vign = inst.slit_vignette()
        assert_quantity_allclose(expected_vign, vign)

    def test_vignette_fwhm2_slit1(self):
        inst_options = {
            "inst_type": "spectrograph",
            "fwhm": 2.0,
            "slit_width": 1.0,
        }

        inst = Instrument(**inst_options)

        expected_vign = 0.434

        vign = inst.slit_vignette(inst_options["slit_width"])

        assert_quantity_allclose(expected_vign, vign)

    def test_vignette_fwhm1_slit1(self):
        inst_options = {
            "inst_type": "spectrograph",
            "fwhm": 1.0,
            "slit_width": 1.0,
        }
        inst = Instrument(**inst_options)

        expected_vign = 0.763

        vign = inst.slit_vignette(inst_options["slit_width"])

        assert_quantity_allclose(expected_vign, vign)

    def test_vignette_fwhm1_slit2(self):
        inst_options = {
            "inst_type": "spectrograph",
            "fwhm": 1.0,
            "slit_width": 2.00,
        }
        inst = Instrument(**inst_options)

        expected_vign = 0.9733

        vign = inst.slit_vignette(inst_options["slit_width"])

        assert_quantity_allclose(expected_vign, vign)

    def test_vignette_fwhm1pt6_slit2(self):
        inst_options = {"inst_type": "spectrograph", "fwhm": 1.6 * u.arcsec, "slit_width": 2.00 * u.arcsec}
        inst = Instrument(**inst_options)

        expected_vign = 0.86125

        vign = inst.slit_vignette(inst_options["slit_width"])

        assert_quantity_allclose(expected_vign, vign)

    def test_vignette_fwhm1_slit2point31(self):
        inst_options = {"inst_type": "spectrograph", "fwhm": 1.0, "slit_width": 2.31}
        inst = Instrument(**inst_options)

        expected_vign = 1.0

        vign = inst.slit_vignette(inst_options["slit_width"])

        assert_quantity_allclose(expected_vign, vign)

    def test_vignette_fwhm2pt0_slit1pt2(self):
        inst_options = {"inst_type": "spectrograph", "fwhm": 2.0 * u.arcsec, "slit_width": 1.20 * u.arcsec}
        inst = Instrument(**inst_options)

        expected_vign = 0.5208

        vign = inst.slit_vignette(inst_options["slit_width"])

        assert_quantity_allclose(expected_vign, vign)

    def test_vignette_fwhm1pt4_fiber2pt0(self):
        inst_options = {
            "inst_type": "spectrograph",
            "fwhm": 1.4 * u.arcsec,
            "fiber_diameter": 2.0 * u.arcsec,
            "slit_width": 2.0 * u.arcsec,
        }
        inst = Instrument(**inst_options)

        expected_vign = 0.7569738014945818

        vign = inst.slit_vignette(inst_options["slit_width"])

        assert_quantity_allclose(expected_vign, vign)

    def test_vignette_fwhm1pt4_fiber2pt0_slit1pt0(self):
        inst_options = {
            "inst_type": "spectrograph",
            "fwhm": 1.4 * u.arcsec,
            "fiber_diameter": 2.0 * u.arcsec,
            "slit_width": 1.0 * u.arcsec,
        }
        inst = Instrument(**inst_options)

        expected_vign = 0.7569738014945818

        vign = inst.slit_vignette(inst_options["slit_width"])

        assert_quantity_allclose(expected_vign, vign)

    def test_vignette_fwhm0pt9_fiber2pt0(self):
        inst_options = {
            "inst_type": "spectrograph",
            "fwhm": 0.9 * u.arcsec,
            "fiber_diameter": 2.0 * u.arcsec,
        }
        inst = Instrument(**inst_options)

        expected_vign = 0.9673838889561505

        vign = inst.slit_vignette()

        assert_quantity_allclose(expected_vign, vign)

    def test_vignette_fwhm2pt0_fiber2pt0(self):
        inst_options = {
            "inst_type": "spectrograph",
            "fwhm": 2.0 * u.arcsec,
            "fiber_diameter": 2.0 * u.arcsec,
        }
        inst = Instrument(**inst_options)

        expected_vign = 0.5

        vign = inst.slit_vignette()

        assert_quantity_allclose(expected_vign, vign)

    def test_vignette_fwhm2pt0_fiber2pt0_nounits(self):
        inst_options = {
            "inst_type": "spectrograph",
            "fwhm": 2.0 * u.arcsec,
            "fiber_diameter": 2.0,
        }
        inst = Instrument(**inst_options)

        expected_vign = 0.5

        vign = inst.slit_vignette()

        assert_quantity_allclose(expected_vign, vign)


class TestMultiChannelInstrument:
    @classmethod
    def setup_class(cls):
        cls.inst_options = {
            "name": "MuSCAT",
            "inst_type": "Imager",
            "ccd_readnoise": 10.0,
            "ccd_gain": 1.0,
            "ccd_xpixels": 2048,
            "ccd_ypixels": 2048,
            "ccd_pixsize": 13.5,
            "ccd_qe": 0.75,
            "channels": {
                "channel1": {"filterlist": ["gp"], "ccd_qe": 0.8},
                "channel2": {"filterlist": ["rp"], "ccd_qe": 0.9},
                "channel3": {"filterlist": ["ip"], "ccd_qe": 0.4},
                "channel4": {"filterlist": ["zs"], "ccd_qe": 0.1, "ccd_gain": 2.0, "ccd_readnoise": 15.0},
            },
        }
        cls.test_mirror_fp = files("tests.etc.data").joinpath("test_mirror.dat").as_posix()

    def test_multichannel_channel_summary(self):
        inst = Instrument(**self.inst_options)

        assert inst.num_channels == 4

    def test_singlechannel_channel_summary(self):
        inst_options = {
            "name": "MuSCAT",
            "inst_type": "Imager",
            "ccd_readnoise": 10.0,
            "ccd_gain": 1.0,
            "ccd_pixsize": 13.5,
        }

        inst = Instrument(**inst_options)

        assert inst.num_channels == 1

    def test_multichannel_ccd_qe(self):
        inst = Instrument(**self.inst_options)

        expected_qe = [0.8, 0.9, 0.4, 0.1]

        assert_quantity_allclose(expected_qe, inst.ccd_qe)

    def test_multichannel_ccd_gain(self):
        inst = Instrument(**self.inst_options)

        expected_gain = [1.0, 1.0, 1.0, 2.0] * u.electron / u.adu

        assert_quantity_allclose(expected_gain, inst.ccd_gain)

    def test_multichannel_ccd_readnoise(self):
        inst = Instrument(**self.inst_options)

        expected_readnoise = [10.0, 10.0, 10.0, 15.0] * u.electron / u.pix

        assert_quantity_allclose(expected_readnoise, inst.ccd_readnoise)

    def test_multichannel_ccd_pixsize(self):
        inst = Instrument(**self.inst_options)

        expected_pixsize = [13.5 * u.micron] * 4

        assert_quantity_allclose(expected_pixsize, inst.ccd_pixsize)

    def test_numchannels(self):
        expected_num_channels = 4

        inst = Instrument(**self.inst_options)

        assert expected_num_channels == inst.num_channels

    def test_modify_numchannels(self):
        expected_num_channels = 4

        inst = Instrument(**self.inst_options)

        assert expected_num_channels == inst.num_channels
        with pytest.raises(AttributeError) as execinfo:
            inst.num_channels = 2
        assert expected_num_channels == inst.num_channels
        assert inst.num_channels != 2

    def test_channel_trans(self):
        inst_options = self.inst_options.copy()
        # Change channel2 to have incorrect ('inst') names
        inst_options["channels"]["channel2"]["num_inst_lenses"] = 1
        inst_options["channels"]["channel2"]["num_inst_mirrors"] = 1
        inst_options["channels"]["channel3"]["trans_components"] = self.test_mirror_fp
        inst_options["channels"]["channel4"]["num_chan_lenses"] = 1
        inst_options["channels"]["channel4"]["num_chan_mirrors"] = 1

        inst = Instrument(**inst_options)

        assert_quantity_allclose(inst.channelset["channel2"].transmission.tpeak(), 1.0)
        assert_quantity_allclose(inst.channelset["channel3"].transmission.tpeak(), 0.92)
        assert_quantity_allclose(inst.channelset["channel4"].transmission.tpeak(), 0.923025)

    def test_channel2filter_map_muscat(self):
        expected_dict = {"channel1": ["gp"], "channel2": ["rp"], "channel3": ["ip"], "channel4": ["zs"]}

        inst = Instrument(**self.inst_options)

        assert expected_dict == inst.channel2filter_map

    def test_channel2filter_map_2filters_per_arm(self):
        expected_dict = {"blue_channel": ["up", "gp"], "red_channel": ["ip", "zs"]}

        inst_options = {
            "name": "Dual Beam Imager",
            "inst_type": "Imager",
            "ccd_readnoise": 10.0,
            "ccd_gain": 1.0,
            "ccd_xpixels": 2048,
            "ccd_ypixels": 2048,
            "ccd_pixsize": 13.5,
            "ccd_qe": 0.75,
            "channels": {
                "blue_channel": {"filterlist": ["up", "gp"], "ccd_qe": 0.5},
                "red_channel": {"filterlist": ["ip", "zs"], "ccd_qe": 0.7},
            },
        }

        inst = Instrument(**inst_options)

        assert expected_dict == inst.channel2filter_map
