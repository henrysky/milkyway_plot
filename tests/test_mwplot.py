import os
import unittest

import numpy.testing as npt

import numpy as np
from astropy import units as  u
import astropy.coordinates as apycoords
from mw_plot import MWSkyMap
from galpy.orbit import Orbit


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


if __name__ == '__main__':
    unittest.main()
