import pytest
import os

from astropy import units as u

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
     assert '<ETC [Site(BPL)]->[Telescope(BPL 1-m)]>' == repr(test_etc)

     assert test_etc.site.name == "BPL"
     assert test_etc.telescope.name == "BPL 1-m"
     assert test_etc.instrument == None

def test_properties():
    test_etc = etc.ETC()

    assert test_etc.site.name == "Maui, HI"
    assert test_etc.site.altitude == 3065 * u.m
    assert test_etc.telescope.name == "FTN"
    assert test_etc.instrument.name == "FLOYDS"
