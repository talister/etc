import pytest
import os
import warnings

import toml
from astropy import units as u
warnings.simplefilter("ignore", pytest.PytestUnknownMarkWarning)
from astropy.tests.helper import assert_quantity_allclose
from synphot.spectrum import SpectralElement

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

    def test_throughput(self):

        optics_options = { 'filterlist' : ['r',] }

        inst = Instrument(**optics_options)

        assert inst.filterlist == ['r', ]
        assert len(inst.filterset) == 1
        assert isinstance(inst.filterset['r'], SpectralElement)
        assert isinstance(inst.throughput('r'), SpectralElement)
        assert_quantity_allclose(inst.throughput('r').tpeak(), u.Quantity(0.78725846), 1e-5)

    def test_throughput_new_ccd_qe(self):

        optics_options = { 'filterlist' : ['r',], 'ccd'  : 0.45 }

        inst = Instrument(**optics_options)

        assert inst.filterlist == ['r', ]
        assert len(inst.filterset) == 1
        assert isinstance(inst.filterset['r'], SpectralElement)
        assert isinstance(inst.throughput('r'), SpectralElement)
        assert_quantity_allclose(inst.throughput('r').tpeak(), u.Quantity(0.78725846/2.0), 1e-5)
