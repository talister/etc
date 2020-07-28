import pytest
from etc import etc


def test_initialize(capsys):
    test_etc = etc.ETC()

    assert ['site', 'telescope', 'instrument'] == test_etc.components
