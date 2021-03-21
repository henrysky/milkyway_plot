import os
import unittest

import numpy.testing as npt

import numpy as np
from astropy import units as  u
import astropy.coordinates as apycoords
from mw_plot import MWSkyMap, MWPlot
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014


class UtilitiesTestCase(unittest.TestCase):
    def test_mw_skymap(self):
        # setup a MWSkyMap instance with projection, other projection can be 'hammer', 'mollweide' etc
        # grayscale: whether to turn the background image to grayscale
        plot_instance = MWSkyMap(projection='aitoff', grayscale=False)

        # so that the colorbar will has a better contract
        # plot_instance.clim = (5., 15.)

        # alpha value for the milkyway image
        plot_instance.imalpha = 1.

        # set up plot title
        plot_instance.title = 'LMC and SMC in red dots'
        plot_instance.s = 200

        # LMC and SMC coordinates, get coordinates with galpy from_name
        lsmc_ra = [Orbit.from_name('LMC').ra(), Orbit.from_name('SMC').ra()] * u.degree
        lsmc_dec = [Orbit.from_name('LMC').dec(), Orbit.from_name('SMC').dec()] * u.degree

        # use mw_scatter instead of scatter
        plot_instance.mw_scatter(lsmc_ra, lsmc_dec, 'r')

        plot_instance.savefig(file='lmc_smc_projection.png')
        
    def test_mw_plot(self):
        # Orbit Integration using galpy for the Sun
        op = Orbit([0., 0., 0., 0., 0., 0.], radec=True, ro=8., vo=220.)
        ts = np.linspace(0, 5, 10000) * u.Gyr
        op.integrate(ts, MWPotential2014)
        x = op.x(ts) * u.kpc
        y = op.y(ts) * u.kpc
        z = op.z(ts)

        # setup a MWPlot instance
        plot_instance = MWPlot(radius=20 * u.kpc, unit=u.kpc, coord='galactocentric', annotation=True)
        plot_instance.imalpha = 1.0
        plot_instance.s = 10  # make the scatter points bigger

        # set up plot title
        plot_instance.title = 'Orbit of Sun in 5Gyr using galpy colored by kpc above galactic plane'

        # use mw_scatter instead of scatter because we want a colorbar
        plot_instance.mw_scatter(x, y, [z, 'kpc above galactic plane'])
        
        plot_instance.savefig(file='gaia.png')


if __name__ == '__main__':
    unittest.main()
