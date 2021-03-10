import pytest
import os
import warnings

import toml
from astropy import units as u
warnings.simplefilter("ignore", pytest.PytestUnknownMarkWarning)
from astropy.tests.helper import assert_quantity_allclose
from synphot.spectrum import SpectralElement, BaseUnitlessSpectrum

from etc.models import *

class TestSite:

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
        test_config_file = toml.load(os.path.abspath(os.path.join(__package__, 'etc', "tests", "data", "test1.toml")))
        site = Site(**test_config_file['site'])

        assert site.name == "BPL"
        assert site.altitude == 17 * u.m
        assert site.latitude == 34.5 * u.deg
        assert site.longitude == -119.86 * u.deg
        assert site.tpeak() == 0.9

    def test_transmission_file(self):
        test_config = { 'name' : "BPL",
                        'altitude' : 7,
                        'latitude' : 34.5,
                        'longitude' : -119.86,
                        'transmission' : os.path.abspath(os.path.join(__package__, 'etc', "tests", "data", "test_atmos.dat"))
                      }
        site = Site(**test_config)

        assert site.name == "BPL"
        assert site.altitude == 7 * u.m
        assert site.latitude == 34.5 * u.deg
        assert site.longitude == -119.86 * u.deg
        assert site.tpeak() == 0.975

    def test_rebin_transmission(self):
        test_config = { 'name' : "BPL",
                        'altitude' : 7,
                        'latitude' : 34.5,
                        'longitude' : -119.86,
                        'transmission' : os.path.abspath(os.path.join(__package__, 'etc', "tests", "data", "test_atmos.dat"))
                      }
        site = Site(**test_config)
        new_waves = u.Quantity([300,600,900,1200], u.nm)

        new_trans = site.rebin_transmission(new_waves)
        assert isinstance(new_trans, SpectralElement)

    def test_sky_brightness(self):
        expected_sky_mags = {'Days_from_New_Moon' : 0, 'U' : 22.0, 'B' : 22.7, 'V' : 21.8, 'R' : 20.9, 'I' : 19.9}

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
        test_config_file = toml.load(os.path.abspath(os.path.join(__package__, 'etc', "tests", "data", "test1.toml")))
        tel = Telescope(**test_config_file['telescope'])

        assert tel.name == "BPL 1-m"
        assert tel.size == 1 * u.m
        assert tel.area == 0.625 * u.m * u.m
        assert tel.num_mirrors == 2
        assert tel.tpeak() == 0.8 * 0.8

    def test_reflectivity_file(self):
        test_config = { 'name' : "BPL 1-m",
                        'size' : 1,
                        'area' : 0.625,
                        'reflectivity' : os.path.abspath(os.path.join(__package__, 'etc', "tests", "data", "test_mirror.dat"))
                      }
        tel = Telescope(**test_config)

        assert tel.name == "BPL 1-m"
        assert tel.size == 1 * u.m
        assert tel.area == 0.625 * u.m * u.m
        assert tel.num_mirrors == 2
        assert tel.tpeak() == 0.92**2

    def test_reflectivity_list(self):
        test_fp = os.path.abspath(os.path.join(__package__, 'etc', "tests", "data", "test_mirror.dat"))

        test_config = { 'name' : "BPL 1-m",
                        'size' : 1,
                        'area' : 0.625,
                        'reflectivity' : [test_fp, test_fp]
                      }
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

    def test_trans_default(self):
        inst = Instrument()

        assert inst.num_mirrors == 0
        assert inst.num_lenses == 1
        assert inst.num_ar_coatings == 2
        assert inst.filterlist == []

        assert inst.transmission.tpeak() == 0.911493 # 0.93 (lens) * 0.99^2 (AR)

    def test_trans_modify_optics(self):
        optics_options = { 'inst_lens_trans' : 0.85,
                           'inst_ar_coating_refl' : 0.95
                         }
        inst = Instrument(**optics_options)


        assert inst.num_mirrors == 0
        assert inst.num_lenses == 1
        assert inst.num_ar_coatings == 2

        assert_quantity_allclose(inst.transmission.tpeak(), 0.767125, 1e-5)

    def test_trans_spectrograph(self):

        optics_options = {  'num_ar_coatings' : 6,    # Air-glass interfaces: prism (2 sides), field flattener (4 sides)
                            'num_inst_lenses' : 3,    # Fused silica prism (two passes) and CCD window
                            'num_inst_mirrors': 4     # Mirrors:  Collimator, Fold, Camera, Grating
                         }
        inst = Instrument(name="FLOYDS", inst_type="SPECTROGRAPH", **optics_options)

        assert inst.name == "FLOYDS"
        assert inst.inst_type == "SPECTROGRAPH"
        assert inst.num_mirrors == 4
        assert inst.num_lenses == 3
        assert inst.num_ar_coatings == 6

        assert_quantity_allclose(inst.transmission.tpeak(), 0.73482187, 1e-5)

    def test_filterset(self):

        optics_options = {  'filterlist' : ['u', 'g', 'r', 'i', 'z'],
                         }
        inst = Instrument(**optics_options)

        assert inst.filterlist == ['u', 'g', 'r', 'i', 'z']
        assert len(inst.filterset) == 5
        assert isinstance(inst.filterset['g'], SpectralElement)

    def test_filterset_primenotation(self):

        optics_options = {  'filterlist' : ['up', 'gp', 'rp', 'ip', 'zs'],
                         }
        inst = Instrument(**optics_options)

        assert inst.filterlist == ['up', 'gp', 'rp', 'ip', 'zs']
        assert len(inst.filterset) == 5
        assert isinstance(inst.filterset['gp'], SpectralElement)
        assert isinstance(inst.filterset['zs'], SpectralElement)

    def test_filterset_mixedcase(self):

        optics_options = { 'filterlist' : ['r', 'R'],
                         }

        inst = Instrument(**optics_options)

        assert inst.filterlist == ['r', 'R']
        assert len(inst.filterset) == 2
        assert isinstance(inst.filterset['r'], SpectralElement)
        assert isinstance(inst.filterset['R'], SpectralElement)
        assert inst.filterset['r'].meta['header']['descrip'] != inst.filterset['R'].meta['header']['descrip']
        assert inst.filterset['r'].meta['expr'] != inst.filterset['R'].meta['expr']

    def test_throughput(self):

        optics_options = { 'filterlist' : ['r',] }

        inst = Instrument(**optics_options)

        assert inst.filterlist == ['r', ]
        assert len(inst.filterset) == 1
        assert isinstance(inst.filterset['r'], SpectralElement)
        assert isinstance(inst.throughput('r'), SpectralElement)
        # 0.93 (lens) * 0.99^2 (AR) * 0.99165802 (filter) * 0.9 (CCD)
        assert_quantity_allclose(inst.throughput('r').tpeak(), 0.81350041, 1e-5)

    def test_throughput_new_ccd_qe(self):

        optics_options = { 'filterlist' : ['r',], 'ccd_qe'  : 0.45 }

        inst = Instrument(**optics_options)

        assert inst.filterlist == ['r', ]
        assert len(inst.filterset) == 1
        assert isinstance(inst.filterset['r'], SpectralElement)
        assert isinstance(inst.throughput('r'), SpectralElement)
        assert_quantity_allclose(inst.throughput('r').tpeak(), 0.81350041/2.0, 1e-5)

    def test_throughput_ccd_qe_quantity(self):

        optics_options = { 'filterlist' : ['r',], 'ccd_qe'  : u.Quantity(0.45, u.dimensionless_unscaled) }

        inst = Instrument(**optics_options)

        assert inst.filterlist == ['r', ]
        assert len(inst.filterset) == 1
        assert isinstance(inst.filterset['r'], SpectralElement)
        assert isinstance(inst.throughput('r'), SpectralElement)
        assert_quantity_allclose(inst.throughput('r').tpeak(), 0.81350041/2.0, 1e-5)

    def test_throughput_ccd_qe_file(self):

        optics_options = { 'filterlist' : ['r',],
                           'ccd_qe'  : os.path.abspath(os.path.join(__package__, 'etc', "tests", "data", "test_ccd_qe.dat"))
                         }

        inst = Instrument(**optics_options)

        assert inst.filterlist == ['r', ]
        assert len(inst.filterset) == 1
        assert isinstance(inst.filterset['r'], SpectralElement)
        assert isinstance(inst.throughput('r'), SpectralElement)
        assert isinstance(inst.ccd_qe, BaseUnitlessSpectrum)
        assert_quantity_allclose(inst.throughput('r').tpeak(), 0.8856121333333334, 1e-5)

    def test_ccd_parameters(self):
        expected_gain = 0.9 * (u.electron / u.adu)
        expected_noise = 4.0 * (u.electron / u.pix)

        ccd_options = { 'ccd_readnoise' : 4.0,
                        'ccd_gain' : 0.9,
                        'ccd_qe'  : os.path.abspath(os.path.join(__package__, 'etc', "tests", "data", "test_ccd_qe.dat"))
                      }

        inst = Instrument(**ccd_options)

        assert isinstance(inst.ccd_qe, BaseUnitlessSpectrum)
        assert_quantity_allclose(inst.ccd_gain, expected_gain)
        assert_quantity_allclose(inst.ccd_readnoise, expected_noise)

    def test_ccdpixsize_default(self):
        expected_value = 17.5 * u.micron

        inst_options = { 'ccd_pixsize' : 17.5 }

        inst = Instrument(**inst_options)

        assert isinstance(inst.ccd_pixsize, u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_pixsize)

    def test_ccdpixsize_goodunits(self):
        expected_value = 13.5 * u.micron

        inst_options = { 'ccd_pixsize' : 0.0135,
                         'ccd_pixsize_units' : "mm"}

        inst = Instrument(**inst_options)

        assert isinstance(inst.ccd_pixsize, u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_pixsize)

    def test_ccdpixsize_goodunits2(self):
        expected_value = 13.5 * u.micron

        inst_options = { 'ccd_pixsize' : 13.5e-6,
                         'ccd_pixsize_units' : "m"}

        inst = Instrument(**inst_options)

        assert isinstance(inst.ccd_pixsize, u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_pixsize)

    def test_ccdpixsize_badunits(self):
        expected_value = 13.5 * u.micron

        inst_options = { 'ccd_pixsize' : 13.5,
                         'ccd_pixsize_units' : "wibble"}

        inst = Instrument(**inst_options)

        assert isinstance(inst.ccd_pixsize, u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_pixsize)


    def test_focalscale_default(self):
        expected_value = 17.5 * (u.arcsec/u.mm)

        inst_options = { 'focal_scale' : 17.5 }

        inst = Instrument(**inst_options)

        assert isinstance(inst.focal_scale, u.Quantity)
        assert_quantity_allclose(expected_value, inst.focal_scale)

    def test_focalscale_goodunits(self):
        expected_value = 1.75 * (u.arcsec/u.cm)

        inst_options = { 'focal_scale' : 1.75,
                         'focal_scale_units' : "arcsec/cm"}

        inst = Instrument(**inst_options)

        assert isinstance(inst.focal_scale, u.Quantity)
        assert_quantity_allclose(expected_value, inst.focal_scale)

    def test_focalscale_badunits(self):
        expected_value = 17.5 * (u.arcsec/u.mm)

        inst_options = { 'focal_scale' : 17.5,
                         'focal_scale_units' : "wibble/cm"}

        inst = Instrument(**inst_options)

        assert isinstance(inst.focal_scale, u.Quantity)
        assert_quantity_allclose(expected_value, inst.focal_scale)

    def test_ccdpixscale(self):
        expected_value = 0.23625 * u.arcsec

        inst_options = { 'focal_scale' : 17.5,
                         'ccd_pixsize' : 13.5}

        inst = Instrument(**inst_options)

        assert isinstance(inst.focal_scale, u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_pixscale)

    def test_ccdpixscale_goodunits(self):
        expected_value = 0.23625 * u.arcsec

        inst_options = { 'focal_scale' : 175,
                         'focal_scale_units' : "arcsec/cm",
                         'ccd_pixsize' : 0.0135,
                         'ccd_pixsize_units' : "mm"}

        inst = Instrument(**inst_options)

        assert isinstance(inst.focal_scale, u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_pixscale)

    def test_ccdpixscale_badunits(self):
        expected_value = 0.23625 * u.arcsec

        inst_options = {
                         'ccd_pixsize' : 13.5,
                         'ccd_pixsize_units' : "wibble",
                         'focal_scale' : 17.5,
                         'focal_scale_units' : "more/wibble"
                         }

        inst = Instrument(**inst_options)

        assert isinstance(inst.focal_scale, u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_pixscale)

    def test_ccdfov_square_defaults(self):
        expected_value = (483.84*u.arcsec, 483.84*u.arcsec)

        inst_options =  {
                          'ccd_pixsize' : 13.5,
                          'focal_scale' : 17.5,
                          'ccd_xpixels' : 2048,
                          'ccd_ypixels' : 2048,
                        }

        inst = Instrument(**inst_options)

        assert isinstance(inst.ccd_fov()[0], u.Quantity)
        assert isinstance(inst.ccd_fov()[1], u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_fov())

    def test_ccdfov_square_arcmin(self):
        expected_value = (483.84/60.*u.arcmin, 483.84/60.0*u.arcmin)

        inst_options =  {
                          'ccd_pixsize' : 13.5,
                          'focal_scale' : 17.5,
                          'ccd_xpixels' : 2048,
                          'ccd_ypixels' : 2048,
                        }

        inst = Instrument(**inst_options)

        fov = inst.ccd_fov(u.arcmin)
        assert isinstance(fov[0], u.Quantity)
        assert isinstance(fov[1], u.Quantity)
        assert_quantity_allclose(expected_value, fov)

    def test_ccdfov_square_badunits(self):
        expected_value = (483.84/60.*u.arcmin, 483.84/60.0*u.arcmin)

        inst_options =  {
                          'ccd_pixsize' : 13.5,
                          'focal_scale' : 17.5,
                          'ccd_xpixels' : 2048,
                          'ccd_ypixels' : 2048,
                        }

        inst = Instrument(**inst_options)

        with pytest.raises(Exception) as execinfo:
            fov = inst.ccd_fov(u.m)

    def test_ccdfov_square_mm(self):
        expected_value = (483.84*u.arcsec, 483.84*u.arcsec)

        inst_options =  {
                          'ccd_pixsize' : 0.0135,
                          'ccd_pixsize_units' : "mm",
                          'focal_scale' : 17.5,
                          'ccd_xpixels' : 2048,
                          'ccd_ypixels' : 2048,
                        }

        inst = Instrument(**inst_options)

        assert isinstance(inst.ccd_fov()[0], u.Quantity)
        assert isinstance(inst.ccd_fov()[1], u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_fov())

    def test_ccdfov_square_mm_cm(self):
        expected_value = (483.84*u.arcsec, 483.84*u.arcsec)

        inst_options =  {
                          'ccd_pixsize' : 0.0135,
                          'ccd_pixsize_units' : "mm",
                          'focal_scale' : 175,
                          'focal_scale_units' : "arcsec/cm",
                          'ccd_xpixels' : 2048,
                          'ccd_ypixels' : 2048,
                        }

        inst = Instrument(**inst_options)

        assert isinstance(inst.ccd_fov()[0], u.Quantity)
        assert isinstance(inst.ccd_fov()[1], u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_fov())

    def test_ccdfov_rect_defaults(self):
        expected_value = (1754*u.arcsec, 1169.3*u.arcsec)

        inst_options =  {
                          'ccd_pixsize' : 9,
                          'ccd_pixsize_units' : "um",
                          'focal_scale' : 63.44,
                          'focal_scale_units' : "arcsec/mm",
                          'ccd_xpixels' : 3072,
                          'ccd_ypixels' : 2048,
                        }

        inst = Instrument(**inst_options)

        assert isinstance(inst.ccd_fov()[0], u.Quantity)
        assert isinstance(inst.ccd_fov()[1], u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_fov(), rtol=1e-1)

    def test_ccdfov_rect_arcmin(self):
        expected_value = (29.2*u.arcmin, 19.5*u.arcmin)

        inst_options =  {
                          'ccd_pixsize' : 9,
                          'ccd_pixsize_units' : "um",
                          'focal_scale' : 63.44,
                          'focal_scale_units' : "arcsec/mm",
                          'ccd_xpixels' : 3072,
                          'ccd_ypixels' : 2048,
                        }

        inst = Instrument(**inst_options)

        assert isinstance(inst.ccd_fov()[0], u.Quantity)
        assert isinstance(inst.ccd_fov()[1], u.Quantity)
        assert_quantity_allclose(expected_value, inst.ccd_fov(u.arcmin), rtol=1e-1)

    def test_vignette_imager(self):
        expected_vign = 1.0

        inst_options = {
                         'inst_type' : 'imager',

                       }

        inst = Instrument(**inst_options)

        vign = inst.slit_vignette()
        assert_quantity_allclose(expected_vign, vign)

    def test_vignette_fwhm2_slit1(self):

        inst_options = {
                         'inst_type' : 'spectrograph',
                         'fwhm' : 2.0,
                         'slit_width' : 1.0,
                       }

        inst = Instrument(**inst_options)

        expected_vign = 0.434

        vign = inst.slit_vignette(inst_options['slit_width'])

        assert_quantity_allclose(expected_vign, vign)

    def test_vignette_fwhm1_slit1(self):

        inst_options = {
                         'inst_type' : 'spectrograph',
                         'fwhm' : 1.0,
                         'slit_width' : 1.0,
                       }
        inst = Instrument(**inst_options)

        expected_vign = 0.763

        vign = inst.slit_vignette(inst_options['slit_width'])

        assert_quantity_allclose(expected_vign, vign)

    def test_vignette_fwhm1_slit2(self):

        inst_options = {
                         'inst_type' : 'spectrograph',
                         'fwhm' : 1.0,
                         'slit_width' : 2.00,
                       }
        inst = Instrument(**inst_options)

        expected_vign = 0.9733

        vign = inst.slit_vignette(inst_options['slit_width'])

        assert_quantity_allclose(expected_vign, vign)

    def test_vignette_fwhm1pt6_slit2(self):

        inst_options = {
                         'inst_type' : 'spectrograph',
                         'fwhm' : 1.6 * u.arcsec,
                         'slit_width' : 2.00 * u.arcsec
                       }
        inst = Instrument(**inst_options)

        expected_vign = 0.86125

        vign = inst.slit_vignette(inst_options['slit_width'])

        assert_quantity_allclose(expected_vign, vign)

    def test_vignette_fwhm1_slit2point31(self):

        inst_options = {
                         'inst_type' : 'spectrograph',
                         'fwhm' : 1.0,
                         'slit_width' : 2.31
                       }
        inst = Instrument(**inst_options)

        expected_vign = 1.0

        vign = inst.slit_vignette(inst_options['slit_width'])

        assert_quantity_allclose(expected_vign, vign)

    def test_vignette_fwhm2pt0_slit1pt2(self):

        inst_options = {
                         'inst_type' : 'spectrograph',
                         'fwhm' : 2.0 * u.arcsec,
                         'slit_width' : 1.20 * u.arcsec
                       }
        inst = Instrument(**inst_options)

        expected_vign = 0.5208

        vign = inst.slit_vignette(inst_options['slit_width'])

        assert_quantity_allclose(expected_vign, vign)
