import pytest
from astroquery.simbad import Simbad


@pytest.fixture(scope="session")
def simbad():
    simbad = Simbad()
    simbad.add_votable_fields("ra(d)", "dec(d)")
    return simbad
