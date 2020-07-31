import pytest
import os

import toml
from astropy import units as u

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

    def test_initialize2(self):
        test_config_file = toml.load(os.path.abspath(os.path.join(__package__, 'etc', "tests", "data", "test1.toml")))
        tel = Telescope(**test_config_file['telescope'])

        assert tel.name == "BPL 1-m"
        assert tel.size == 1 * u.m
        assert tel.area == 0.625 * u.m * u.m
        assert tel.num_mirrors == 2
