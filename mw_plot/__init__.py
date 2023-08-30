try:
    import bokeh
    _BOKEH_LOADED = True
except ImportError:
    _BOKEH_LOADED = False
else:
    from mw_plot.mw_plot_bokeh import *

from mw_plot.mw_plot_matplotlib import *
from pkg_resources import get_distribution
from astropy.coordinates import SkyCoord, ICRS
import astropy.units as u
import numpy as np

version = __version__ = get_distribution("mw_plot").version

# RA and DEC of galactic center, galactic anti-center, galactic north and south pole in degree
center_radec = [266.4167, -29.0078]
anti_center_radec = [86.4167, 28.0078]
northpole_radec = [192.7667, 27.1167]
southpole_radec = [12.7667, -27.1167]


def mw_radec(deg=True, size=3600):
    """
    Get RA DEC coordinates of the milkyway

    :param deg: To get RA/DEC in deg when True or rad when False
    :type deg: bool
    :param size: number of point
    :type size: int

    :return: RA, DEC
    :rtype: np.ndarray

    :History: 2021-Feb-26 - Written - Henry Leung (University of Toronto)
    """
    c = SkyCoord(
        np.linspace(0.0, 360.0, size) * u.deg, np.zeros(size) * u.deg, frame="galactic"
    )
    c = c.transform_to(ICRS)
    idx = np.argsort(c.ra.to(u.deg).value)
    if deg:
        return c.ra.to(u.deg).value[idx], c.dec.to(u.deg).value[idx]
    else:
        return c.ra.to(u.rad).value[idx], c.dec.to(u.rad).value[idx]
