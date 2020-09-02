import pytest
import os
import warnings

import toml
from astropy import units as u
warnings.simplefilter("ignore", pytest.PytestUnknownMarkWarning)
from astropy.tests.helper import assert_quantity_allclose
from synphot.spectrum import SpectralElement, BaseUnitlessSpectrum

from etc.models import *

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
        assert tel.tpeak() == 0.91**3

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

        assert inst.transmission.tpeak() == 0.88209

    def test_trans_modify_optics(self):
        optics_options = { 'inst_lens_trans' : 0.85,
                           'inst_ar_coating_refl' : 0.95
                         }
        inst = Instrument(**optics_options)


        assert inst.num_mirrors == 0
        assert inst.num_lenses == 1
        assert inst.num_ar_coatings == 2

        assert inst.transmission.tpeak() == 0.767125

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

        assert inst.transmission.tpeak() == 0.6659793414426959

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
        assert_quantity_allclose(inst.throughput('r').tpeak(), 0.78725846, 1e-5)

    def test_throughput_new_ccd_qe(self):

        optics_options = { 'filterlist' : ['r',], 'ccd_qe'  : 0.45 }

        inst = Instrument(**optics_options)

        assert inst.filterlist == ['r', ]
        assert len(inst.filterset) == 1
        assert isinstance(inst.filterset['r'], SpectralElement)
        assert isinstance(inst.throughput('r'), SpectralElement)
        assert_quantity_allclose(inst.throughput('r').tpeak(), 0.78725846/2.0, 1e-5)

    def test_throughput_ccd_qe_quantity(self):

        optics_options = { 'filterlist' : ['r',], 'ccd_qe'  : u.Quantity(0.45, u.dimensionless_unscaled) }

        inst = Instrument(**optics_options)

        assert inst.filterlist == ['r', ]
        assert len(inst.filterset) == 1
        assert isinstance(inst.filterset['r'], SpectralElement)
        assert isinstance(inst.throughput('r'), SpectralElement)
        assert_quantity_allclose(inst.throughput('r').tpeak(), 0.78725846/2.0, 1e-5)

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
        assert_quantity_allclose(inst.throughput('r').tpeak(), 0.857044, 1e-5)

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
