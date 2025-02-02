import pytest
from astroquery.simbad import Simbad


@pytest.fixture(scope="session")
def simbad():
    simbad = Simbad()
    return simbad
