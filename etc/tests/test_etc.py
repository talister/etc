import pytest
import os
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
