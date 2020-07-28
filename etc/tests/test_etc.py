import pytest
from etc import etc


def test_initialize():
    test_etc = etc.ETC()

    assert ['site', 'telescope', 'instrument'] == test_etc.components

def test_initialize_specified_components():
    test_etc = etc.ETC(components=['sky', 'earth', 'underground cavern'])

    assert ['sky', 'earth', 'underground cavern'] == test_etc.components
