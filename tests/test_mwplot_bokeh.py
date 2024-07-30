import numpy as np
from astropy import units as u
from mw_plot import MWSkyMapBokeh, MWPlotBokeh
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014


def test_mw_skymap_bokeh():
    # setup a MWSkyMapBokeh instance
    # grayscale: whether to turn the background image to grayscale
    plot_instance = MWSkyMapBokeh(grayscale=False)

    # so that the colorbar will has a better contract
    # plot_instance.clim = (5., 15.)

    # alpha value for the milkyway image
    plot_instance.imalpha = 1.0

    # set up plot title
    plot_instance.title = "LMC and SMC in red dots"
    plot_instance.s = 200

    # LMC and SMC coordinates, get coordinates with galpy from_name
    lsmc_ra = [Orbit.from_name("LMC").ra(), Orbit.from_name("SMC").ra()] * u.degree
    lsmc_dec = [
        Orbit.from_name("LMC").dec(),
        Orbit.from_name("SMC").dec(),
    ] * u.degree

    # use mw_scatter instead of scatter
    plot_instance.scatter(lsmc_ra, lsmc_dec)

def test_mw_plot_bokeh():
    # Orbit Integration using galpy for the Sun
    op = Orbit([0.0, 0.0, 0.0, 0.0, 0.0, 0.0], radec=True, ro=8.0, vo=220.0)
    ts = np.linspace(0, 5, 10000) * u.Gyr
    op.integrate(ts, MWPotential2014)
    x = op.x(ts) * u.kpc
    y = op.y(ts) * u.kpc
    z = op.z(ts)

    # setup a MWPlotBokeh instance
    plot_instance = MWPlotBokeh(
        radius=20 * u.kpc, unit=u.kpc, coord="galactocentric", annotation=True
    )
    plot_instance.imalpha = 1.0
    plot_instance.s = 10  # make the scatter points bigger

    # set up plot title
    plot_instance.title = (
        "Orbit of Sun in 5Gyr using galpy colored by kpc above galactic plane"
    )

    # use mw_scatter instead of scatter because we want a colorbar
    plot_instance.scatter(x, y)
